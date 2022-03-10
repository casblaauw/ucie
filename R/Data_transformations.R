Sys.setenv("DISPLAY"=":0.0")

#' @import dplyr
#' @import rgl
#' @import remotes
#' @import colorspace
#' @importFrom geometry convhulln
#' @importFrom geometry delaunayn
#' @importFrom pracma distmat
#' @importFrom ptinpoly pip3d
NULL

# RGL_USE_NULL=TRUE


RGB_space <- as.matrix(data.frame(
  "R"= c(seq(0, 255, by=32),255, # K -> R
         seq(0, 255, by=32),255, # G -> Y
         rep(255, 9), # Y -> W
         rep(255, 9), # R -> Y
         rep(255, 9), # R -> M
         seq(0, 255, by=32),255, # B -> M
         rep(0, 9), # B -> C
         seq(0, 255, by=32),255, # C -> W
         rep(0, 9), # G -> C
         rep(0, 9), # K -> G
         rep(0, 9), # K -> B
         rep(255, 9)), # M -> W

  "G"= c(rep(0, 9), # K -> R
         rep(255, 9), # G -> Y
         rep(255, 9), # Y -> W
         seq(0, 255, by=32),255, # R -> Y
         rep(0, 9),  # R -> M
         rep(0, 9),  # B -> M
         seq(0, 255, by=32),255, # B -> C
         rep(255, 9), # C -> W
         rep(255, 9), # G -> C
         seq(0, 255, by=32),255, # K -> G
         rep(0, 9), # K -> B
         seq(0, 255, by=32),255), # M -> W

  "B"= c(rep(0, 9), # K -> R
         rep(0, 9), # G -> Y
         seq(0, 255, by=32),255, # Y -> W
         rep(0, 9), # R -> Y
         seq(0, 255, by=32),255, # R -> M
         rep(255, 9),  # B -> M
         rep(255, 9), # B -> C
         rep(255, 9), # C -> W
         seq(0, 255, by=32),255, # G -> C
         rep(0, 9), # K -> G
         seq(0, 255, by=32),255, # K -> B
         rep(255, 9))# M -> W
))

#' Check and reformat data for use with U-CIE functions
#'
#' Takes a 2D or 3D matrix or data frame and returns a 3D matrix, checked for non-numeric values.
#'
#' @param dataset A data frame or matrix of 2D or 3D values.
#' @param check_3D Optional: whether to check and force 3D output. Needed for reduced/cielab data,
#' not needed for KNN data.
#'
#' @details
#'   The function does a number of processing steps to accept as many data formats as feasible:
#'   * Rownames in the first column of the matrix (as tested by unconvertability to numeric) or of the data frame are discarded.
#'   * Character matrices are converted to numeric.
#'   * Data frames are converted to matrices.
#'   * If rownames_col is given, that column is discarded as well.
#'   * 2D matrices or data frames are expanded to 3D by padding with 1's.
#'
#' @keywords internal
PrepData <- function(dataset, check_3D = TRUE) {

  # If data frame, transform to matrix
  if (inherits(dataset, "data.frame")) {
    # If first column looks like rownames, move to actual rownames
    if (is.character(dataset[[1]])) {
      rownames(dataset) <- dataset[[1]]
      dataset <- dataset[[-1]]
    }

    # Check for non-numeric data
    if (!all(sapply(dataset, is.numeric))) {
      stop(paste("The dataset contains non-numeric columns.",
                         "Please check your data and use the rownames_col argument if necessary."))
    }
    dataset <- as.matrix(dataset)
  }

  # If rownames encoded in first matrix row, move to true rownames if not present or drop
  if (anyNA(as.numeric(dataset[,1]))) {
    val <- dataset[,1]
    dataset <- dataset[,-1]
    if (is.null(rownames(dataset))) {
      rownames(dataset) <- val
    }
  }

  # If character matrix, transform into numeric
  if (!all(is.numeric(dataset))) {
    class(dataset) <- 'numeric'
  }

  # Check for missing/leftover character data
  if (anyNA(dataset)) {
    stop(paste("Your matrix data contains missing data or",
                       "non-numeric data outside the first column.",
                       "Please check your dataset."))
  }

  if (check_3D) {

    # Check data dimensionality
    if (ncol(dataset) == 2) {
      warning("Data expanded to 3D!")
      dataset <- cbind(dataset, 1)
    }
    if (ncol(dataset) != 3) {
      stop("The dataset should have 3 numeric columns!")
    }

    # Check for missing values
    if (anyNA(dataset)) {
      stop("The dataset has missing values. Check again!")
    }

    # Set column names
    colnames(dataset) <- c('X', 'Y', 'Z')
  }

  return(dataset)
}

#' Convert RGB coordinates to LAB coordinates
#'
#' @param data tibble, data frame or matrix with columns `R`/`G`/`B`,
#' scaled 0-1 or 0-255
#'
#' @return A data frame or tibble (if supplied), otherwise a matrix, with columns
#' `L`/`A`/`B`.
#'
#' @keywords internal
RGB2Lab <- function(data) {
  if (inherits(data, 'tbl_df')) {
    output_class <- 'tibble'
  } else if (inherits(data, 'data.frame')) {
    output_class <- 'data.frame'
  } else {
    output_class <- 'matrix'
  }
  data <- as.matrix(data)
  if (any(data > 1)) {
    data <- data / 255
  }
  RGB_obj <- colorspace::RGB(R = data[,'R'], G = data[,'G'], B = data[,'B'])
  LAB_obj <- as(RGB_obj, 'LAB')
  LAB_coords <- colorspace::coords(LAB_obj)
  if (output_class == 'tibble') {
    return(dplyr::as_tibble(LAB_coords))
  } else if (output_class == 'data.frame') {
    return(as.data.frame(LAB_coords))
  } else {
    return(LAB_coords)
  }
}

#' Return CIELAB outline space coordinates and colours
#'
#' @param RGB_points RGB data to move to UCIE. Default is set of points
#' on RGB extremes cube edges with step size 32.
#'
#' @return A data frame with columns L/A/B/colour, containing coordinates and hex codes.
#' @keywords internal
CIELABspace <- function(RGB_points = RGB_space){
  CIELAB_coords <- RGB2Lab(RGB_points)
  CIELAB <- dplyr::bind_cols(
    CIELAB_coords,
    colour = colorspace::hex(colorspace::LAB(CIELAB_coords), fixup = TRUE)
  )
  return(CIELAB)
}

DataConvex <- function(Query){
  ch_cloud <- geometry::convhulln(as.matrix(Query))
  ConvexCloud <- as.matrix(Query)[unique(as.vector(ch_cloud)),] # Only vertices of convex hull
  return(ConvexCloud)
}

Rotation <- function(ConvexCloud, RotL, Rota, Rotb){
  if (is.na(RotL) | is.na(Rota) | is.na(Rotb)) {
    print(paste('An rotation value is NA, which will likely lead to an error. Values:', RotL, Rota, Rotb))
  }
  ConvexCloud <-  suppressWarnings(rgl::rotate3d(obj = ConvexCloud, angle = RotL, x = 1, y = 0, z = 0))
  ConvexCloud <-  suppressWarnings(rgl::rotate3d(obj = ConvexCloud, angle = Rota, x = 0, y = 1, z = 0))
  ConvexCloud <-  suppressWarnings(rgl::rotate3d(obj = ConvexCloud, angle = Rotb, x = 0, y = 0, z = 1))
  return(ConvexCloud)
}

Translation <- function(ConvexCloud, TrL, Tra, Trb){
  ConvexCloud <-  suppressWarnings(rgl::translate3d(ConvexCloud, TrL, Tra, Trb))
  return(ConvexCloud)
}

Scaling <- function(ConvexCloud, S){
  ConvexCloud <-  suppressWarnings(rgl::scale3d(ConvexCloud, S, S, S))
  return(ConvexCloud)
}

TransformedConvexCloud <- function(S, RotL, Rota, Rotb,  TrL, Tra, Trb, WL, Wa, Wb, Query){
  ConvexCloud <- Rotation(Query, RotL, Rota, Rotb)
  ConvexCloud <- Scaling(ConvexCloud, S)
  ConvexCloud <- Translation(ConvexCloud, TrL, Tra, Trb)

  L <- (max(ConvexCloud[,1]) - min(ConvexCloud[,1]))
  a <- (max(ConvexCloud[,2]) - min(ConvexCloud[,2]))
  b <- (max(ConvexCloud[,3]) - min(ConvexCloud[,3]))

  if (b > L & b > a) {
    ConvexCloud <- ConvexCloud
  } else {
    ConvexCloud[,1] <- WL*ConvexCloud[,1]
    ConvexCloud[,2] <- Wa*ConvexCloud[,2]
    ConvexCloud[,3] <- Wb*ConvexCloud[,3]
  }
  return(ConvexCloud)
}

#' Calculate distance of outside points to colour space.
#'
#' Determines which points are outside the colour space convex hull, calculates
#' the distance matrix between those points and the colour space, then, for each
#' data point, returns the distance to the closest colour space hull point.
#'
#' @param S Scaling parameter.
#' @param RotL L* rotation parameter.
#' @param Rota a* rotation parameter.
#' @param Rotb b* rotation parameter.
#' @param TrL L* translation parameter.
#' @param Tra a* translation parameter.
#' @param Trb b* translation parameter.
#' @param WL Weight of L* axis.
#' @param Wa Weight of a* axis.
#' @param Wb Weight of b* axis.
#' @param Query (Convex hull of the) data to fit, as a matrix.
#' @param polygon Vertices (points) of the colour space convex hull, as a matrix.
#' Like output returned by [DataConvex()].
#' @param faces Edges (lines) of the colour space convex hull, as a matrix.
#' Like output returned by [geometry::convhulln()].
#'
#' @return A vector with distance to the colour space, for each point outside the space,
#' or 0 if no points are outside the space.
#' @keywords internal

Distance <- function(S, RotL, Rota, Rotb,  TrL, Tra, Trb, WL, Wa, Wb, Query, polygon, faces){
  ConvexCloud <- TransformedConvexCloud(S, RotL, Rota, Rotb,  TrL, Tra, Trb, WL, Wa, Wb, Query)
  point_in_space <- ptinpoly::pip3d(polygon, faces, ConvexCloud)
  outside <- ConvexCloud[which(point_in_space==-1), ]
  if(length(outside) == 0){
    dist <- 0
  } else {
    dist_mat <- pracma::distmat(polygon, outside)
    dist <- as.data.frame(apply(dist_mat,2,min))
  }
  return(dist)
} # Gives me the distance of the points from the Polygon

#' Data fit minimalisation function
#'
#' Takes fit and weighting parameters, the convex hull of the data, and the
#' convex hull (vertices and faces) of the colour space to fit the data to.
#' Returns a fit metric, consisting of the scaling parameter, penalised strongly
#' by any points outside the colour space.
#'
#' The penalty function calculates the distance of all outside points to the
#' closest colour space hull point, then subtracts the sum of all squared distances.
#'
#' @param param Vector of fitting parameters (S, RotL, Rota, Rotb, TrL, Tra, Trb)
#' @param WL Weight of L* axis.
#' @param Wa Weight of a* axis.
#' @param Wb Weight of b* axis.
#' @param data (Convex hull of the) data to fit, as a matrix.
#' @param polygon Vertices (points) of the colour space convex hull, as a matrix.
#' Like output returned by [DataConvex()].
#' @param faces Edges (lines) of the colour space convex hull, as a matrix.
#' Like output returned by [geometry::convhulln()].
#'
#' @keywords internal
ObjectiveFunction <- function(param, WL, Wa, Wb, data, polygon, faces) {
  outside_dist <- Distance(param[1], param[2], param[3], param[4], param[5], param[6], param[7], WL, Wa, Wb, data, polygon, faces)
  f <- param[1] - sum(outside_dist^2)
  return(f)
}

TranslationGuess <- function(dataset, colorspace) {
  centroidval_colorspace <- PolygonCentroid(colorspace) # centroid of color space
  centroidval_cloud <- PolygonCentroid(dataset) # centroid of cloud

  guess <- centroidval_colorspace - centroidval_cloud %>%
    setNames(c('TrL', "Tra", "Trb"))

  return(guess)
}

ScalingGuess <- function(dataset, polygon) {
  #--- Scaling factor ---#
  ConvexCloud <- (DataConvex(dataset))

  Cloud1Size <- max(ConvexCloud[, 1]) - min(ConvexCloud[, 1])
  Cloud2Size <- max(ConvexCloud[, 2]) - min(ConvexCloud[, 2])
  Cloud3Size <- max(ConvexCloud[, 3]) - min(ConvexCloud[, 3])

  # Find the smallest and use it as L

  Cloud_scaled <- matrix(nrow = nrow(ConvexCloud), ncol = ncol(ConvexCloud))
  if (Cloud1Size < Cloud2Size & Cloud1Size < Cloud3Size) {
    MaxScalingFactor_1 <- max(polygon[,1]) / Cloud1Size
    MaxScalingFactor_2 <- max(polygon[,2]) / Cloud2Size
    MaxScalingFactor_3 <- max(polygon[,3]) / Cloud3Size

    MaxScalingFactor <- ifelse(MaxScalingFactor_1 < MaxScalingFactor_2,MaxScalingFactor_1, MaxScalingFactor_2)
    MaxScalingFactor <- ifelse(MaxScalingFactor < MaxScalingFactor_3, MaxScalingFactor, MaxScalingFactor_3)

    # #offset
    # Cloud1offset <- (max(ConvexCloud[, 1]) + min(ConvexCloud[, 1])) / 2
    # Cloud2offset <- (max(ConvexCloud[, 2]) + min(ConvexCloud[, 2])) / 2
    # Cloud3offset <- (max(ConvexCloud[, 3]) + min(ConvexCloud[, 3])) / 2
    #
    # Cloud_scaled[, 1] <- (ConvexCloud[, 1] - Cloud1offset) * MaxScalingFactor + 50
    # Cloud_scaled[, 2] <- (ConvexCloud[, 2] - Cloud2offset) * MaxScalingFactor
    # Cloud_scaled[, 3] <- (ConvexCloud[, 3] - Cloud3offset) * MaxScalingFactor
  }
  if (Cloud2Size < Cloud1Size & Cloud2Size < Cloud3Size) {
    MaxScalingFactor_1 <- max(polygon[,2]) / Cloud1Size
    MaxScalingFactor_2 <- max(polygon[,1]) / Cloud2Size
    MaxScalingFactor_3 <- max(polygon[,3]) / Cloud3Size

    MaxScalingFactor <- ifelse(MaxScalingFactor_1 < MaxScalingFactor_2, MaxScalingFactor_1, MaxScalingFactor_2)
    MaxScalingFactor <- ifelse(MaxScalingFactor < MaxScalingFactor_3, MaxScalingFactor, MaxScalingFactor_3)

    # #offset
    # Cloud1offset <- (max(ConvexCloud[, 1]) + min(ConvexCloud[, 1])) / 2
    # Cloud2offset <- (max(ConvexCloud[, 2]) + min(ConvexCloud[, 2])) / 2
    # Cloud3offset <- (max(ConvexCloud[, 3]) + min(ConvexCloud[, 3])) / 2
    #
    # Cloud_scaled[, 1] <- (ConvexCloud[, 2] - Cloud2offset) * MaxScalingFactor + 50
    # Cloud_scaled[, 2] <- (ConvexCloud[, 1] - Cloud1offset) * MaxScalingFactor
    # Cloud_scaled[, 3] <- (ConvexCloud[, 3] - Cloud3offset) * MaxScalingFactor
  }
  if (Cloud3Size < Cloud1Size & Cloud3Size < Cloud2Size) {
    MaxScalingFactor_1 <- max(polygon[,3]) / Cloud1Size
    MaxScalingFactor_2 <- max(polygon[,2]) / Cloud2Size
    MaxScalingFactor_3 <- max(polygon[,1]) / Cloud3Size

    MaxScalingFactor <- ifelse( MaxScalingFactor_1 < MaxScalingFactor_2, MaxScalingFactor_1, MaxScalingFactor_2)
    MaxScalingFactor <- ifelse(MaxScalingFactor < MaxScalingFactor_3, MaxScalingFactor, MaxScalingFactor_3)

    # #offset
    # Cloud1offset <- (max(ConvexCloud[, 1]) + min(ConvexCloud[, 1])) / 2
    # Cloud2offset <- (max(ConvexCloud[, 2]) + min(ConvexCloud[, 2])) / 2
    # Cloud3offset <- (max(ConvexCloud[, 3]) + min(ConvexCloud[, 3])) / 2
    #
    # Cloud_scaled[, 1] <- (ConvexCloud[, 3] - Cloud3offset) * MaxScalingFactor  + 50
    # Cloud_scaled[, 2] <- (ConvexCloud[, 1] - Cloud1offset) * MaxScalingFactor
    # Cloud_scaled[, 3] <- (ConvexCloud[, 2] - Cloud2offset) * MaxScalingFactor
  }

  return(c(S = MaxScalingFactor))
}

#' Calculate centroid of the data's convex hull polygon
#'
#' Calculates the centroid of a given data matrix's convex hull, using Delaunay triangulation.
#'
#' Effectively finds the centroids of the tetrahedrons returned by Delaunay triangulation
#' and their volume. The weighted mean of those centroids, weighted by volume, is
#' the polygon's centroid. Unlike the arithmetic mean of the convex hull's vertices,
#' it is not biased towards more 'curved' sides that have more vertices.
#'
#' @param dataset A 3D matrix of either full data, or only the points making up
#' that data's convex hull (the vertices, like from [DataConvex()]).
#'
#' @keywords internal
PolygonCentroid <- function(dataset) {
  del <- geometry::delaunayn(dataset, full = TRUE)

  # Calculate centroids of all delaunay tetrahedrons:
  # Take mean of x/y/z coord respectively of each of the 4 coordinates per row
  del_centroids <- cbind(
    x = rowMeans(cbind(dataset[del$tri[,1], 1], dataset[del$tri[,2], 1], dataset[del$tri[,3], 1], dataset[del$tri[,4], 1])),
    y = rowMeans(cbind(dataset[del$tri[,1], 2], dataset[del$tri[,2], 2], dataset[del$tri[,3], 2], dataset[del$tri[,4], 2])),
    z = rowMeans(cbind(dataset[del$tri[,1], 3], dataset[del$tri[,2], 3], dataset[del$tri[,3], 3], dataset[del$tri[,4], 3]))
  )
  del_weights <- del$areas

  centroid <- apply(del_centroids, 2, weighted.mean, del_weights)
  return(setNames(centroid, colnames(dataset)))
}

#' Fit data to CIELAB space
#'
#' @param dataset Dataset to fit, as a matrix.
#' @param WL Weight of L* axis.
#' @param Wa Weight of a* axis.
#' @param Wb Weight of b* axis.
#' @param center Boolean: whether to center the data and color space during the
#' fitting process, which improves the fit performance but means the data needs
#' to be centered before the returned parameters can be used. Default is TRUE.
#'
#' @return A named vector of fitting parameters (S, RotL, Rota, Rotb, TrL, Tra, Trb).
#' @keywords internal
FitColorsFunction <- function(dataset, WL, Wa, Wb, center = TRUE){
  LAB_space <- RGB2Lab(RGB_space)

  if (center) {
    data_centroid <- PolygonCentroid(dataset)
    dataset <- Translation(dataset, -data_centroid[1], -data_centroid[2], -data_centroid[3])

    LAB_centroid <- PolygonCentroid(LAB_space)
    LAB_space <- Translation(LAB_space, -LAB_centroid[1], -LAB_centroid[2], -LAB_centroid[3])
  }

  dat_polygon <- DataConvex(dataset)
  LAB_polygon <- DataConvex(LAB_space)

  LAB_chull <- geometry::convhulln(LAB_polygon, output.options = 'FA')


  #------ Initial guess ---------------------------------------------------------#
  S_guess <- ScalingGuess(dat_polygon, LAB_polygon)
  dat_polygon_scaled <- Scaling(dat_polygon, S_guess)

  # Create a list of start values, forcibly mirroring L-rotation to explore rotation landscape
  start_params <- lapply(
    1:25,
    function(counter) {
      Rot <- c(RotL = pi / 4 + (pi / (1 + 0.5 * counter)),
               Rota = pi / 4,
               Rotb = pi / 4)

      Tr <- dat_polygon_scaled %>%
        Rotation(Rot[1], Rot[2], Rot[3]) %>%
        TranslationGuess(LAB_space)

      params <- c('S' = 1, Rot, Tr) #S, RotL, Rota, Rotb, TrL, Tra, Trb

      return(params)
    }
  )
  #------ Simplex optimizer -----------------------------------------------------#

  # set.seed(123)

  # Optimise for each set of start params and return in data frame
  fitted_params <- lapply(start_params,
                           function(params) {
                             optim_result <- stats::optim(
                               par = params,
                               method = "Nelder-Mead",
                               ObjectiveFunction,
                               WL = WL,
                               Wa = Wa,
                               Wb = Wb,
                               data = dat_polygon_scaled,
                               polygon = LAB_polygon,
                               faces = LAB_chull$hull,
                               control = list(fnscale = -1, maxit = 1000)
                             )
                             return(dplyr::tibble(
                               'params' = list(optim_result$par),
                               'value' = optim_result$value
                             ))
                           })
  fitted_params <- dplyr::bind_rows(fitted_params)

  result <- fitted_params$params[[which.max(fitted_params$value)]]
  result[1] <- S_guess * result[1]
  return(result)
}



#' Mapping 3D Data into CIELab color space
#'
#' @param dataset 3-column dataset to be translated into colors.
#' @param WL Weight of L* axis in optimization function. Default value 1.
#' @param Wa Weight of a* axis in optimization function. Default value 1.
#' @param Wb Weight of b* axis in optimization function. Default value 1.
#' @param S Scaling factor for color mapping. Default value 1.
#' @param LAB_coordinates Logical. If FALSE, the function returns a data frame with hex colors.
#' If TRUE, the function returns a data frame with the L*a*b* coordinates. Default value FALSE.
#' @param center Boolean: whether to center the data and color space during the
#' fitting process. Generally improves fit performance. Default is TRUE.
#' @param verbose Boolean: whether to print progress updates. Default is FALSE.
#'
#' @return Returns a data frame with 2 or 4 columns: a rownames column and either
#' hex colour codes or L/A/B coordinates.
#' @export
#'
#' @examples
#' \donttest{
#'   df <- data.frame(V1 = runif(10, 0, 1), V2 = runif(10, 0, 5), V3 = runif(10, 0, 30))
#'   data2cielab(df, Wb = 1.2, S = 1.6)
#'   data2cielab(df, LAB_coordinates = TRUE)
#' }
data2cielab <- function(dataset, WL = 1, Wa = 1, Wb = 1, S = 1, LAB_coordinates = FALSE, center = TRUE, verbose = FALSE){

  dataset <- PrepData(dataset)

  if (!is.null(rownames(dataset))) {
    data_rownames <- rownames(dataset)
  } else {
    data_rownames <- 1:nrow(dataset)
  }

  if (verbose) start_fit <- Sys.time()
  final_params <- FitColorsFunction(dataset, WL, Wa, Wb, center = center)
  if (verbose) {
    end_fit <- Sys.time()
    print('Fitting data to CIELAB space:')
    print(end_fit - start_fit)
    print('Fitted params:')
    print(final_params)
    print('Transforming data')
  }

  # Transform data
  dataset <- Scaling(dataset, final_params[1]*S)
  dataset <- Rotation(dataset, final_params[2], final_params[3], final_params[4])
  dataset <- Translation(dataset, final_params[5], final_params[6], final_params[7])

  if (center) {
    if (verbose) print('Un-centering data')
    # Move data from center to CIELAB space
    data_centroid <- PolygonCentroid(dataset)
    cielab_centroid <- PolygonCentroid(RGB2Lab(RGB_space))
    translation_set <- cielab_centroid - data_centroid
    dataset <- Translation(dataset, translation_set[1], translation_set[2], translation_set[3])
  }

  if (verbose) print('Turning data into colours')
  colorspace_obj <- colorspace::LAB(round(dataset, 2))

  if (LAB_coordinates) {
    # Turn coordinates matrix into a data frame
    col_coords <- colorspace::coords(colorspace_obj) %>%
      as.data.frame()
    # Move row names or indices to column for tidy data
    col_coords <- cbind(names = data_rownames, col_coords)
    rownames(col_coords) <- NULL
    return(col_coords)
  } else {
    # Turn hex character vector into a data frame
    col_hex <- colorspace::hex(colorspace_obj, fixup = TRUE) %>%
      data.frame(colour = .)
    # Move row names or indices to column for tidy data
    col_hex <- cbind(names = data_rownames, col_hex)
    rownames(col_hex) <- NULL
    return(col_hex)
  }
}








