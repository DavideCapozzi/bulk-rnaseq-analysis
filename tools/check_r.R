# Pre-commit check for R files, run from the repository root.
#   style: fails if styler would change a file (nothing is rewritten)
#   lint:  fails on any lintr finding (settings from .lintr)
# Usage: Rscript tools/check_r.R style|lint FILE...
args <- commandArgs(trailingOnly = TRUE)
mode <- args[[1]]
files <- args[-1]
if (length(files) == 0) quit(status = 0)

if (mode == "style") {
  changed <- Filter(function(f) {
    out <- styler::style_file(f, dry = "on", include_roxygen_examples = TRUE)
    isTRUE(out$changed)
  }, files)
  if (length(changed) > 0) {
    message("styler would change: ", paste(changed, collapse = ", "), "\nrun: Rscript -e 'styler::style_pkg()'")
    quit(status = 1)
  }
} else if (mode == "lint") {
  # object_usage_linter resolves the package's internal functions through its namespace: load it from the source,
  # so the result does not depend on whether (or which version of) the package is installed.
  pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
  lints <- lapply(files, lintr::lint)
  n <- sum(lengths(lints))
  if (n > 0) {
    for (l in lints) if (length(l)) print(l)
    quit(status = 1)
  }
} else {
  stop("usage: Rscript tools/check_r.R style|lint FILE...")
}
