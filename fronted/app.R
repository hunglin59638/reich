#!/usr/bin/env Rscript
library(shiny)
cmdArgs <- commandArgs(trailingOnly = FALSE)
root_dir <- normalizePath(dirname((sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)]))))[1]
port <- cmdArgs[length(cmdArgs)]
if (is.na(port) || grepl("^--file=", port, fixed = F)) {
    port <- 59638
} else {
    port <- as.numeric(port)
}
options(shiny.port = port)
options(shiny.host = "0.0.0.0")
runApp(appDir = root_dir, launch.browser = F)
