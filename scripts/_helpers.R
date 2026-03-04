\
# Helper functions used across scripts

safe_pdf <- function(filename, width = 10, height = 5) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  pdf(filename, width = width, height = height)
}

save_session_info <- function(out_file = "results/logs/session_info.txt") {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  si <- capture.output(sessionInfo())
  writeLines(si, con = out_file)
}

set_analysis_seed <- function(seed = 1234) {
  set.seed(seed)
  if ("future" %in% (.packages())) {
    # future uses RNG streams; this ensures consistent behavior when possible
    options(future.rng.onMisuse = "ignore")
  }
}
