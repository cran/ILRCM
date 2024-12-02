#' @title Simulated Irregular Longitudinal Data
#' @name sdata
#' @docType data
#' @description Simulated irregular longitudinal data for 1000 patients.
#' This dataset contains irregularly spaced time points and responses for analysis.
#' @usage data(sdata)
#' @format A data frame with 8631 rows and 3 variables:
#' \describe{
#'   \item{subject_id}{ID of subjects}
#'   \item{time}{Irregular time points.}
#'   \item{response}{Response values at different time points.}
#' }
#' @examples
#' data(sdata)
#' head(sdata)
"sdata"
