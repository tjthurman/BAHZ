#' Make a prior configuration file
#'
#' Makes a configuration file which you can use to specify the priors for your
#' hybrid zone analysis. By default, creates the file in the working directoy
#' with the name \code{prior_config_template.yaml}. You can supply a filepath
#' and a different name, if you want. Will not overrite files by default, must
#' set overwrite = T.
#'
#' DETAILS TO BE ADDED.
#'
#' @param path The folder where you wish to save the prior configuration file.
#'   The default is the current working directory.
#' @param name The name you want to give the prior configuration file. The
#'   default is \code{prior_config_template.yaml}
#' @param overwrite Do you want to overwrite existing files?
#'   \code{TRUE} or \code{FALSE}, default is \code{FALSE}
#'
#' @return The function returns a value of NULL. \code{invidible(NULL)}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Make config file in current working directory with default name.
#' # Won't overwrite
#' make_prior_config()
#'
#' # Make config file in a new path with a new name.
#' # Will overwrite
#' make_prior_config(path = "/path/to/directory", name = "my_priors.yaml", overwrite = T)
#' }
#' # TO ADD
#'

make_prior_config <- function(path = getwd(),
                                  name = "prior_config_template.yaml",
                                  overwrite = F) {
  assertthat::assert_that(length(name) == 1, msg = "name must be of length 1")
  assertthat::assert_that(is.character(name) == T, msg = "name must be a character string")
  assertthat::assert_that(stringr::str_detect(name, pattern = "\\.yaml$") == 1, msg = "name must end in .yaml")

  assertthat::assert_that(is.logical(overwrite) == T, msg = "overwrite must be either TRUE (T) or FALSE (F)")


  newfile <- file.path(normalizePath(path), name, fsep = .Platform$file.sep)
  if (overwrite == F) { # If overwriting is not allowed
    if (file.exists(newfile) == T) { # Stop if file exists
      stop("Template file already exists, set overwrite = T is you want to overwrite it")
    }
    else {
      file.copy(from = system.file("extdata", "prior_config_template.yaml", package = "cap", mustWork = TRUE),
                to = newfile,
                overwrite = overwrite)
    }
  }
  else { # If overwriting is allowed
    if (file.exists(newfile) == T) { # Notify when overwriting
      message(paste0(name, " alrady exists, overwriting"))
      file.copy(from = system.file("extdata", "prior_config_template.yaml", package = "cap", mustWork = TRUE),
                to = newfile,
                overwrite = overwrite)
    } else {
      file.copy(from = system.file("extdata", "prior_config_template.yaml", package = "cap", mustWork = TRUE),
                to = newfile,
                overwrite = overwrite)
    }
  }
  if (file.exists(newfile) == T) {
    message("Configuration file successfully generated")
  } else {
    stop("Configuration file not generated")
  }
  return(invisible(NULL))
}



