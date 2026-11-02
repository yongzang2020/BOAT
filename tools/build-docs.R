#!/usr/bin/env Rscript
# Run from any directory: Rscript /path/to/BOAT/tools/build-docs.R
# Requires msm, knitr, rmarkdown, Pandoc, and a LaTeX installation.
# The original vignette's eval flags and simulation settings are preserved.

build_docs <- function() {
  script_arg <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_arg) != 1L) stop("Run this file with Rscript.")
  script_path <- sub("^--file=", "", script_arg)
  # R encodes spaces in its --file argument as ~+~ on some platforms.
  if (!file.exists(script_path)) script_path <- gsub("~+~", " ", script_path, fixed = TRUE)
  script_file <- normalizePath(script_path, mustWork = TRUE)
  root <- dirname(dirname(script_file))
  dependencies <- c("msm", "knitr", "rmarkdown")
  missing <- dependencies[!vapply(dependencies, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install these R packages first: ", paste(missing, collapse = ", "))

  if (dir.exists("/Library/TeX/texbin")) {
    Sys.setenv(PATH = paste("/Library/TeX/texbin", Sys.getenv("PATH"), sep = .Platform$path.sep))
  }
  if (!nzchar(Sys.which("pdflatex"))) stop("pdflatex is required to build the PDFs.")
  if (!rmarkdown::pandoc_available()) {
    arch <- if (grepl("arm64|aarch64", R.version$arch)) "aarch64" else "x86_64"
    candidates <- c(
      file.path("/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools", arch),
      "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools",
      "/Applications/RStudio.app/Contents/MacOS/pandoc"
    )
    for (candidate in candidates) {
      if (file.exists(file.path(candidate, "pandoc"))) {
        Sys.setenv(RSTUDIO_PANDOC = candidate)
        break
      }
    }
    rmarkdown::find_pandoc(cache = FALSE)
  }
  if (!rmarkdown::pandoc_available()) stop("Pandoc is required to build the guide.")

  work <- tempfile("boat-docs-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  snapshot <- file.path(work, "BOAT")
  lib <- file.path(work, "library")
  output <- file.path(work, "output")
  dir.create(snapshot)
  dir.create(lib)
  dir.create(output)
  package_items <- c("DESCRIPTION", "NAMESPACE", "R", "man", "vignettes", "inst", "data", "src", "LICENSE", "LICENCE")
  package_items <- package_items[file.exists(file.path(root, package_items))]
  if (!all(file.copy(file.path(root, package_items), snapshot, recursive = TRUE))) {
    stop("Could not copy package files into the temporary build directory.")
  }
  old_lib <- .libPaths()
  on.exit(.libPaths(old_lib), add = TRUE)
  .libPaths(c(lib, old_lib))
  old_wd <- setwd(work)
  on.exit(setwd(old_wd), add = TRUE)
  rbin <- file.path(R.home("bin"), "R")
  run_r_command <- function(arguments, log_name) {
    result <- system2(rbin, arguments, stdout = TRUE, stderr = TRUE)
    writeLines(result, file.path(work, log_name))
    status <- attr(result, "status")
    if (!is.null(status) && status != 0L) {
      stop(paste(tail(result, 30L), collapse = "\n"))
    }
  }
  run_r_command(c("CMD", "INSTALL", "--no-multiarch", "--no-test-load",
                  paste0("--library=", shQuote(lib)), shQuote(snapshot)), "install.log")
  stopifnot(normalizePath(find.package("BOAT")) == normalizePath(file.path(lib, "BOAT")))
  run_r_command(c("CMD", "Rd2pdf", "--no-preview", "--force",
                  paste0("--output=", shQuote(file.path(output, "BOAT-manual.pdf"))),
                  shQuote(snapshot)), "manual.log")

  options(width = 80)
  knitr::opts_chunk$set(size = "small")
  header <- file.path(work, "guide-header.tex")
  writeLines(c(
    "\\usepackage{fvextra}",
    "\\widowpenalty=10000",
    "\\clubpenalty=10000",
    "\\DefineVerbatimEnvironment{Highlighting}{Verbatim}{breaklines,commandchars=\\\\\\{\\}}"
  ), header)
  format <- rmarkdown::pdf_document(
    toc = TRUE, number_sections = FALSE, latex_engine = "pdflatex",
    includes = rmarkdown::includes(in_header = header)
  )
  rmarkdown::render(
    file.path(snapshot, "vignettes", "BOAT.Rmd"), output_format = format,
    output_file = "BOAT-guide.pdf", output_dir = output,
    intermediates_dir = file.path(work, "render"),
    envir = new.env(parent = globalenv()), quiet = FALSE
  )
  docs <- file.path(root, "docs")
  dir.create(docs, showWarnings = FALSE)
  files <- c("BOAT-manual.pdf", "BOAT-guide.pdf")
  if (!all(file.copy(file.path(output, files), docs, overwrite = TRUE))) {
    stop("Could not copy the completed PDFs to docs/.")
  }
  message("Created: ", paste(file.path(docs, files), collapse = "\n"))
}

build_docs()
