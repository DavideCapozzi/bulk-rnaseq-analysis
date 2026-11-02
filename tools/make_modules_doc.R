# Writes docs/MODULES.md from inst/modules.yaml and inst/schema/schema.yaml, run from the repository root.
# Usage: Rscript tools/make_modules_doc.R
pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
writeLines(bulkrnaseq:::modules_doc(), "docs/MODULES.md")
