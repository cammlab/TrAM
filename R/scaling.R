#'
#' \code{median_abs_diff} returns the median absolute difference between adjacent data values
#'
#' @param values A numeric vector of time ordered data values grouped by \code{index}
#' @param index A vector of grouping indexes defining distinct time series instances
#' @return A numeric value
#' @examples
#'   num.vecs <- 10
#'   vec.lengths <- sample(10:15, 10, replace=TRUE)
#'   index <- rep(1:num.vecs, vec.lengths)
#'   values <- rep(100*rnorm(num.vecs), vec.lengths) + rnorm(sum(vec.lengths))
#'   median_abs_diff(values, index)
#' 
median_abs_diff <- function(values, index) {
    stopifnot(all(!is.na(values))) # may not have any NA
    stopifnot(is.numeric(values)) # must be numeric
    stopifnot(length(values) > 1) # must have more than one value
    
    data_list <- split(values, index) # split up the data by index, keeping the order within each vector

    ## for each individual vector, compute the absolute value of each time difference
    all_diffs <- unlist(lapply(data_list, function(timeseries) abs(diff(timeseries))))

    ## compute and return the grand median of all those time differences
    median(all_diffs)
}

#'
#' Return a data frame which contains time series from an input data frame rescaled
#' by each time series's median absolute temporal fluctuations
#'
#' @param data A data frame containing a time series grouping column, a time
#'   column (assumed uniformly spaced), and multiple time series columns
#' @param group.col.name Name of the column in \code{data} which contains the
#'   time series grouping
#' @param time.col.name Name of the column in \code{data} which contains the
#'   time for each row
#' @param data.col.names Character vector of column names containing simultaneous
#'   time series values
#' @return A data frame containing the same column names specified in the function
#'   call, but with the data column values rescaled
#' @examples
#'   num.time.series = 10
#'   time.series.length = 25
#'   num.data.points <-  num.time.series * time.series.length
#'   data = data.frame(group=rep(1:num.time.series, time.series.length),
#'                     time=rep(1:time.series.length, num.time.series),
#'                     ts1=rnorm(num.data.points),
#'                     ts2=10*rnorm(num.data.points))
#'   rescaled.data <-  median_abs_diff_rescale(data, "group", "time", c("ts1", "ts2"))
#' 
median_abs_diff_rescale <- function(data, group.col.name, time.col.name, data.col.names) {
    ## arrange the data by group and then by time so the data are ordered
    data <- data[order(data[,group.col.name], data[,time.col.name]),]

    index <- data[,group.col.name] # data vector ids
    time <- data[,time.col.name]
    
    rescaled.values <- lapply(data.col.names, function(name) {
        values <- data[,name]
        scaling.factor <- median_abs_diff(values, index)

        values/scaling.factor # rescale
    })

    result.frame <-  as.data.frame(cbind(index, time, do.call(cbind, rescaled.values)))
    names(result.frame) <-  c(group.col.name, time.col.name, data.col.names)

    result.frame
}

