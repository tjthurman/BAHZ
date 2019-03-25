.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  for (m in modules) loadModule(m, what = TRUE)
}

.onAttach <- function(...) {
  packageStartupMessage("Running bahz in parallel may speed up model fitting.\n",
                        "See startup message from Rstan, below:\n",
                        "\n",
                        "For execution on a local, multicore CPU with excess RAM we recommend calling\n",
                        "options(mc.cores = parallel::detectCores())")
}
