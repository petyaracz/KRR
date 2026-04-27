#' krr: Kernel Ridge Regression for Phonological Analogy
#'
#' Fit and apply kernel ridge regression models using precomputed phonological
#' distances. Key functions: [krr::train_krr()], [krr::predict_krr()], [krr::check_krr_inputs()].
#'
#' See `vignette("lakok", package = "krr")` for a worked example.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom dplyr all_of arrange filter mutate select slice
#' @importFrom tidyr crossing pivot_wider
#' @importFrom purrr map2_dbl
#' @importFrom tibble tibble
#' @importFrom stats cor plogis qlogis
NULL
