#!/usr/bin/env Rscript

# Compatibility entry point for older local instructions.  The maintained
# implementation is main.R; sourcing it forwards the same command-line
# arguments and therefore enforces the three-input workflow contract.
script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_file_arg[[1]]), mustWork = FALSE))
} else {
  getwd()
}
source(file.path(script_dir, "main.R"))
