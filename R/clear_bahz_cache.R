#' Clear cache of saved results from memoised functions in bahz
#'
#' Clear the cache of all already calculated results for all the memoised functions
#' (predict_cline, cline_summary) in bahz.
#'
#' @import memoise
#'
#' @return invisible(NULL)
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' clear_bahz_cache()
#'
#' }
#'
clear_bahz_cache <- function() {
  memoise::forget(bahz::predict_cline)
  memoise::forget(bahz::cline_summary)
  print("bahz cache cleared")
  invisible(NULL)
}
