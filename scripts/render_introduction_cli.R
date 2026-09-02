#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
proj_root <- if (length(args)) args[[1]] else getwd()

source(file.path(proj_root, "scripts", "render_introduction.R"), local = TRUE)
render_introduction(proj_root = proj_root)
