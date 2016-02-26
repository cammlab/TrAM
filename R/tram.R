#' Compute the Tracking Aberration Measure (TrAM)
#'
#' \code{tram} returns the Tracking Aberration Measure computed from multiple time series
#'
#' This function provides measures of discontinuities present across multiple simultaneous
#' time series to enable detection of erroneous data that are "jumpy".
#'
#' @param values A numeric matrix whose columns are simultaneous time series. Must
#'   have column names and at least 6 time points (rows). Values from different
#'   types of measurement will typically be rescaled so their typical fluctuations
#'   are comparable in scale (see \code{\link{median_abs_diff_rescale}}).
#' @param num.knots A positive integer number of knots in the smoothing spline fits.
#' @param p Positive exponent used to combine jumpiness across the time series.
#'   Values >1 tend to emphasize the single time series with the most jumpiness,
#'   whereas values <1 emphasize simultaneous jumps between multiple time series.
#' @param euclidian.list A named list whose elements are character vectors. Each
#'   character vector specifies a group of column names whose deviations from
#'   smoothness will be combined in a Euclidian (exponent=2) manner prior to
#'   combining with other time series. This is useful for positional coordinates
#'   where spherical symmetry in the error is important to maintain.
#' @return Vector containing overall TrAM (named "tram") and component TrAM
#'   values prior to summing.
#' @examples
#'
#' num.objects <- 100
#' num.timepoints <- 25
#' jumpy.fraction <- 0.1
#' jump.scale <- 1
#' noise.level <- 0.2
#'
#' data <- do.call(rbind, lapply(1:num.objects, function(object.number) {
#'     generate_smooth_series <- function() {
#'         times.scaled <- 2*(1:num.timepoints)/num.timepoints - 1
#'         rnorm(1) + rnorm(1)*times.scaled + 3*rnorm(1)*times.scaled^2 # quadratic
#'     }
#'
#'     relative.brightness.scale <- 100
#'
#'     x <- generate_smooth_series() + noise.level*rnorm(num.timepoints)
#'     y <- generate_smooth_series() + noise.level*rnorm(num.timepoints)
#'     brightness <- relative.brightness.scale * (10 + generate_smooth_series() +
#'                                                    noise.level*rnorm(num.timepoints))
#'
#'     ## see if this one will be jumpy
#'     is.jumpy <- runif(1) < jumpy.fraction
#'     if(is.jumpy) {
#'         jump.timepoint <- 2 + round((num.timepoints-5)*runif(1))
#'
#'         x[jump.timepoint] <-  x[jump.timepoint] + sample(c(-1,1),1) * jump.scale*(1+0.3*rnorm(1))
#'         y[jump.timepoint] <-  y[jump.timepoint] + sample(c(-1,1),1) * jump.scale*(1+0.3*rnorm(1))
#'         brightness[jump.timepoint] <-  brightness[jump.timepoint] +
#'             relative.brightness.scale*sample(c(-1,1),1)*jump.scale*(1+0.3*rnorm(1))
#'     }
#'
#'     data.frame(object.number, time=1:num.timepoints, x, y, brightness, is.jumpy)
#' }))
#'
#' scaled.data <- median_abs_diff_rescale(data, "object.number", "time", c("x", "y", "brightness"))
#' trams <- do.call(rbind, by(scaled.data, scaled.data$object.number, function(scaled.object.data) {
#'     values <- with(scaled.object.data, cbind(x, y, brightness))
#'     data.frame(object.number=scaled.object.data$object.number[1], tram=tram(values)['tram'])
#' }))
#'
#' trams <- merge(trams, subset(data, time == 1, # merge in is.jumpy annotation
#'  select=c(object.number, is.jumpy)))
#'
#'
#' \dontrun{
#' with(trams, plot(tram, is.jumpy)) # dot plot
#' }
#'
#' @export
tram <- function(values,
                 num.knots = floor(nrow(values)/5),
                 p = 0.5,
                 euclidian.list = list()) {

    if(nrow(values) < 6) stop("Must be at least 6 points in time series.")
    stopifnot(!is.null(colnames(values)))
    stopifnot(is.numeric(values))
    stopifnot(num.knots > 0)
    stopifnot(p > 0)


    ## if there are any NA then we can't do the calculation
    if(any(is.na(values))) return(NA)

    time.points <- 1:nrow(values) # define uniformly spaced time points for the values

    ## get a smoothed version of each trajectory
    splined.values <- apply( t(values), 1, function(v) smooth.spline(time.points, v, nknots=num.knots)$y)

    abs.deltas <- abs(splined.values - values) # absolute differences from smoothed
    num.time.series <- ncol(abs.deltas) # number of dimensions


    ## If there are no Euclidian metrics then that makes things easy
    if(length(euclidian.list) == 0) {
        weights <- rep (1/num.time.series, num.time.series) # all variables have the same weight
        names(weights) <- colnames(values)
    } else {
        cols.list <- list() # where to save up columns
        used.names <- vector("character")

        weights <- vector("numeric") # to store weights that are equal to the number of initial variables in each final variable

        ## loop through each euclidian vector by name
        euclidian.var.names <- names(euclidian.list)
        for(euclidian.var.name in euclidian.var.names) {
            v <- euclidian.list[[euclidian.var.name]] # get the names of the data that form this Euclidian metric
            sm <- abs.deltas[, v] # get those data values
            result.col <- as.matrix(sqrt(rowMeans(sm^2))) # now compute the Euclidian "distance"
            colnames(result.col) <- euclidian.var.name # rename it
            cols.list[[length(cols.list)+1]] <- result.col

            weights[euclidian.var.name] <- length(v) # set the weight to be the number of initial variables
            used.names <- c(used.names, v) # keep track of all names already used
        }

        remaining.col.names <- setdiff(colnames(abs.deltas), used.names) # names not used for euclidian metrics

        if(length(remaining.col.names) > 0) { # if there are any left over take account of them
            remaining.abs.deltas <- abs.deltas[, remaining.col.names, drop=F] # need to keep as matrix so use drop=F
            cols.list[[length(cols.list)+1]] <- remaining.abs.deltas # store it at the end
            weights[colnames(remaining.abs.deltas)] <- 1 # all remaining columns have weight=1
        }

        weights <- weights / sum(weights) # normalize the weights

        ## now get the replacement matrix
        abs.deltas <- do.call(cbind, cols.list)
    }

    v <- (rowSums(weights[colnames(abs.deltas)] * abs.deltas^p)) ^ (1/p)
    multi.tram <- max(v) # dist exponent < 1 means conincidences across dimensions weigh more heavily

    individual.trams <- apply(abs.deltas, 2, max) # individual tram is max for that column (across time)

    result <- c(multi.tram, individual.trams)
    names(result) <- paste0(names(result), ".tram") # give each datum a nice name
    names(result)[1] <- "tram" # base name for the multi tram

    return(result)
}

