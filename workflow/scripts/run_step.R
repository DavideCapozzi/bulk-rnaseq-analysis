#!/usr/bin/env Rscript
# Thin wrapper: Rscript run_step.R <step> <config> <repo> [module]
# Loads the package from the repository checkout, so the workflow always runs the code of the commit recorded
# in the run manifest, then delegates to bulkrnaseq::run_step_cli().
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: run_step.R <step> <config> <repo> [module]")
suppressMessages(pkgload::load_all(args[[3]], export_all = FALSE, helpers = FALSE, quiet = TRUE))
bulkrnaseq::run_step_cli(args[[1]], args[[2]], module = if (length(args) >= 4) args[[4]] else NULL)
