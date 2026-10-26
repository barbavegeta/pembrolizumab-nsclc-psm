#' Read the model inputs table
#'
#' @param path Path to `data-raw/inputs.csv`.
#' @return The inputs as a data frame. Every row must name a source.
#' @export
read_inputs <- function(path) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (anyDuplicated(x$parameter)) stop("Duplicated parameter names in ", path)
  missing_src <- x$parameter[is.na(x$source) | x$source == ""]
  if (length(missing_src) > 0) {
    stop("Inputs without a source: ", paste(missing_src, collapse = ", "))
  }
  x
}

#' Look up one input value
#'
#' @param inputs Output of [read_inputs()].
#' @param name Parameter name.
#' @param numeric Convert the value to numeric (default `TRUE`).
#' @return The value.
#' @export
input_value <- function(inputs, name, numeric = TRUE) {
  v <- inputs$value[inputs$parameter == name]
  if (length(v) != 1) stop("Input '", name, "' not found exactly once")
  if (numeric) as.numeric(v) else v
}
