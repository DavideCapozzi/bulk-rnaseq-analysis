#!/usr/bin/env Rscript

################################################################################
# R Package Installation Script 
################################################################################

suppressPackageStartupMessages({
  options(warn = 1, error = function() {
    cat("\n[FATAL ERROR] Installation failed. Check r_install_compact.log
")
    q(status = 1)
  })
})

################################################################################
# CONFIGURATION
################################################################################

BIOC_VERSION <- "3.22"
CRAN_MIRROR <- "https://cloud.r-project.org"
LOG_FILE <- "r_install_compact.log"
ERROR_REPORT <- "r_install_errors.txt"

# Package lists with priority levels
SKIP_PACKAGES <- Sys.getenv("SKIP_PACKAGES", "")  # User can set this
VERBOSE_MODE <- Sys.getenv("VERBOSE_R_INSTALL", "false") == "true"

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

################################################################################
# LOGGING (MINIMAL CONSOLE, FULL FILE)
################################################################################

log_file <- function(msg, level = "INFO") {
  ts <- timestamp()
  full_msg <- sprintf("[%s] [%s] %s", ts, level, msg)
  cat(full_msg, "
", file = LOG_FILE, append = TRUE)
  if (VERBOSE_MODE) cat(full_msg, "
")
}

log_console <- function(msg) {
  # Only important messages to console
  cat(msg)
  cat(msg, "
", file = LOG_FILE, append = TRUE)
}

################################################################################
# SETUP PHASE
################################################################################

setup_environment <- function() {
  cat("╔═══════════════════════════════════════════════════════════════════════════╗")
  cat("║         R Package Installation - Compact v3                               ║")
  cat("╚═══════════════════════════════════════════════════════════════════════════╝")

  log_file("Installation started at", "INFO")
  log_file(paste("R Version:", R.version$version.string), "INFO")
  log_file(paste("Working directory:", getwd()), "INFO")

  # Check library paths
  lib_paths <- .libPaths()
  log_file(paste0("Library paths: ", paste(lib_paths, collapse = ", ")), "INFO")

  # Ensure main library is writable
  main_lib <- lib_paths[1]
  if (file.access(main_lib, 2) != 0) {
    log_console("✗ ERROR: Main library path is not writable!")
    log_console(paste0("  Path: ", main_lib))
    q(status = 1)
  }

  log_console("✓ Environment ready")
  log_console("")
}

################################################################################
# BIOCMANAGER SETUP
################################################################################

setup_biocmanager <- function() {
  log_console("[SETUP] Installing/loading BiocManager...")

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    tryCatch({
      install.packages("BiocManager", repos = CRAN_MIRROR, quiet = TRUE)
    }, error = function(e) {
      log_console("✗ Failed to install BiocManager")
      log_file(paste("BiocManager install error:", e$message), "ERROR")
      q(status = 1)
    })
  }

  library(BiocManager, quietly = TRUE, warn.conflicts = FALSE)

  # Set version
  tryCatch({
    options(BioC_mirror = "https://bioconductor.org")
    BiocManager::install(version = BIOC_VERSION, ask = FALSE, update = FALSE)
  }, error = function(e) {
    log_file(paste("BiocManager version setup error:", e$message), "WARNING")
  })

  log_console("✓ BiocManager ready")
  log_console("")
}

################################################################################
# INSTALL PACKAGE WITH RETRY AND ERROR HANDLING
################################################################################

install_pkg <- function(pkg, source = "bioc", retry = FALSE) {
  # Check if should skip
  if (grepl(pkg, SKIP_PACKAGES)) {
    log_file(paste0(pkg, " - SKIPPED (user specified)"), "INFO")
    return(list(success = NA, reason = "User skipped"))
  }

  # Check if already installed
  if (requireNamespace(pkg, quietly = TRUE)) {
    ver <- packageVersion(pkg)
    log_file(paste0(pkg, " - Already installed (v", ver, ")"), "SKIP")
    return(list(success = TRUE, version = as.character(ver)))
  }

  # Attempt installation
  tryCatch({
    if (source == "bioc") {
      BiocManager::install(pkg,
                          version = BIOC_VERSION,
                          ask = FALSE,
                          update = FALSE,
                          quiet = TRUE,
                          force = retry)
    } else {
      install.packages(pkg, repos = CRAN_MIRROR, quiet = TRUE, dependencies = TRUE)
    }
  }, error = function(e) {
    log_file(paste0(pkg, " - Installation error: ", e$message), "ERROR")
    return(NULL)
  }, warning = function(w) {
    log_file(paste0(pkg, " - Warning: ", w$message), "WARNING")
    return(NULL)
  })

  # Verify installation
  if (requireNamespace(pkg, quietly = TRUE)) {
    ver <- packageVersion(pkg)
    log_file(paste0(pkg, " - SUCCESS (v", ver, ")"), "SUCCESS")
    return(list(success = TRUE, version = as.character(ver)))
  } else {
    log_file(paste0(pkg, " - Verification FAILED"), "ERROR")

    # Try to get error details
    tryCatch({
      loadNamespace(pkg)
    }, error = function(e) {
      log_file(paste0(pkg, " - Load error: ", e$message), "ERROR")
    })

    return(list(success = FALSE, reason = "Verification failed"))
  }
}

################################################################################
# INSTALLATION PHASES
################################################################################

install_all_packages <- function() {
  log_console("[PHASE 1] Critical Base Dependencies...")

  critical <- c("BiocGenerics", "S4Vectors", "IRanges", "XVector", 
                "GenomeInfoDb", "zlibbioc")

  critical_results <- list()
  for (pkg in critical) {
    result <- install_pkg(pkg, "bioc")
    critical_results[[pkg]] <- result
    status <- if (result$success) "✓" else "✗"
    cat(sprintf("  %s %s
", status, pkg), file = LOG_FILE, append = TRUE)
  }

  log_console(paste0("  Completed: ", sum(sapply(critical_results, function(x) x$success)), 
                     "/", length(critical_results), " packages"))
  log_console("")

  # Core packages
  log_console("[PHASE 2] Core DE Analysis Packages...")

  core <- list(
    c("tximport", "bioc"), c("DESeq2", "bioc"), c("edgeR", "bioc"),
    c("limma", "bioc"), c("GenomicFeatures", "bioc"),
    c("SummarizedExperiment", "bioc"), c("Biostrings", "bioc")
  )

  core_results <- list()
  for (item in core) {
    result <- install_pkg(item[1], item[2])
    core_results[[item[1]]] <- result
    status <- if (result$success) "✓" else if (is.na(result$success)) "⊘" else "✗"
    cat(sprintf("  %s %s
", status, item[1]), file = LOG_FILE, append = TRUE)
  }

  success_count <- sum(sapply(core_results, function(x) isTRUE(x$success)))
  log_console(paste0("  Completed: ", success_count, "/", length(core_results), " packages"))
  log_console("")

  # Annotation packages
  log_console("[PHASE 3] Annotation Packages...")

  annot <- list(
    c("AnnotationDbi", "bioc"), c("org.Hs.eg.db", "bioc"),
    c("GO.db", "bioc"), c("KEGG.db", "bioc")
  )

  annot_results <- list()
  for (item in annot) {
    result <- install_pkg(item[1], item[2])
    annot_results[[item[1]]] <- result
    status <- if (result$success) "✓" else if (is.na(result$success)) "⊘" else "✗"
    cat(sprintf("  %s %s
", status, item[1]), file = LOG_FILE, append = TRUE)
  }

  success_count <- sum(sapply(annot_results, function(x) isTRUE(x$success)))
  log_console(paste0("  Completed: ", success_count, "/", length(annot_results), " packages"))
  log_console("")

  # Enrichment packages
  log_console("[PHASE 4] Enrichment Analysis Packages...")

  enrich <- list(
    c("fgsea", "bioc"), c("clusterProfiler", "bioc"), c("enrichplot", "bioc"),
    c("DOSE", "bioc"), c("pathview", "bioc")
  )

  enrich_results <- list()
  for (item in enrich) {
    result <- install_pkg(item[1], item[2])
    enrich_results[[item[1]]] <- result
    status <- if (result$success) "✓" else if (is.na(result$success)) "⊘" else "✗"
    cat(sprintf("  %s %s
", status, item[1]), file = LOG_FILE, append = TRUE)
  }

  success_count <- sum(sapply(enrich_results, function(x) isTRUE(x$success)))
  log_console(paste0("  Completed: ", success_count, "/", length(enrich_results), " packages"))
  log_console("")

  # Visualization
  log_console("[PHASE 5] Visualization & Utilities...")

  viz <- list(
    c("ggplot2", "cran"), c("pheatmap", "cran"), c("RColorBrewer", "cran"),
    c("viridis", "cran"), c("scales", "cran"), c("reshape2", "cran"),
    c("cowplot", "cran"), c("ggrepel", "cran"), c("gplots", "cran"),
    c("readr", "cran"), c("readxl", "cran"), c("writexl", "cran"),
    c("msigdbr", "cran")
  )

  viz_results <- list()
  for (item in viz) {
    result <- install_pkg(item[1], item[2])
    viz_results[[item[1]]] <- result
    status <- if (result$success) "✓" else if (is.na(result$success)) "⊘" else "✗"
    cat(sprintf("  %s %s
", status, item[1]), file = LOG_FILE, append = TRUE)
  }

  success_count <- sum(sapply(viz_results, function(x) isTRUE(x$success)))
  log_console(paste0("  Completed: ", success_count, "/", length(viz_results), " packages"))
  log_console("")

  # Reproducibility
  log_console("[PHASE 6] Reproducibility Tools...")

  repro <- list(
    c("renv", "cran"), c("rmarkdown", "cran"), c("knitr", "cran")
  )

  repro_results <- list()
  for (item in repro) {
    result <- install_pkg(item[1], item[2])
    repro_results[[item[1]]] <- result
    status <- if (result$success) "✓" else if (is.na(result$success)) "⊘" else "✗"
    cat(sprintf("  %s %s
", status, item[1]), file = LOG_FILE, append = TRUE)
  }

  success_count <- sum(sapply(repro_results, function(x) isTRUE(x$success)))
  log_console(paste0("  Completed: ", success_count, "/", length(repro_results), " packages"))
  log_console("")

  # Combine all results
  all_results <- c(critical_results, core_results, annot_results, 
                   enrich_results, viz_results, repro_results)

  return(all_results)
}

################################################################################
# SUMMARY REPORT
################################################################################

print_summary <- function(results) {
  total <- length(results)
  success <- sum(sapply(results, function(x) isTRUE(x$success)))
  failed <- sum(sapply(results, function(x) !is.na(x$success) && !x$success))
  skipped <- sum(sapply(results, function(x) is.na(x$success)))

  cat("╔═══════════════════════════════════════════════════════════════════════════╗")
  cat("║                         INSTALLATION SUMMARY                              ║")
  cat("╚═══════════════════════════════════════════════════════════════════════════╝")
  log_console(sprintf("  Total packages:    %3d", total))
  log_console(sprintf("  ✓ Successful:      %3d", success))
  log_console(sprintf("  ✗ Failed:          %3d", failed))
  log_console(sprintf("  ⊘ Skipped:         %3d", skipped))
  log_console(sprintf("  Success Rate:      %.1f%%", 100 * success / (success + failed)))
  log_console("")

  if (failed > 0) {
    log_console("Failed packages:")
    for (pkg_name in names(results)) {
      result <- results[[pkg_name]]
      if (!is.na(result$success) && !result$success) {
        log_console(paste0("  ✗ ", pkg_name, " - ", result$reason))
      }
    }
    log_console("")
    log_console("To retry failed packages:")
    log_console("  SKIP_PACKAGES='' Rscript env/install_r_packages_v3_compact.r")
    log_console("")
    log_console("Or to skip and continue with successful packages:")
    failed_pkg_list <- paste(names(results)[!sapply(results, function(x) isTRUE(x$success))], 
                            collapse = "|")
    log_console(paste0("  SKIP_PACKAGES='", failed_pkg_list, "' Rscript env/install_r_packages_v3_compact.r"))
  }

  log_console("")
  log_console("For detailed information, see:")
  log_console(paste0("  - Log file: ", LOG_FILE))
  log_console("")
}

################################################################################
# MAIN EXECUTION
################################################################################

setup_environment()
setup_biocmanager()

results <- install_all_packages()
print_summary(results)

q(save = "no", status = 0)
