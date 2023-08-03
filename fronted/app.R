#!/usr/bin/env Rscript
library(shiny)
library(yaml)
library(tidyverse)
cmdArgs <- commandArgs(trailingOnly = FALSE)
REICH_ROOT <- normalizePath(dirname(
    (sub("^--file=", "", cmdArgs[grep("^--file=", cmdArgs)]))
))[1] %>% dirname()

CONFIG <- read_yaml(file.path(REICH_ROOT, "config.yml"))

port <- cmdArgs[length(cmdArgs)]
if (is.na(port) || grepl("^--file=", port, fixed = F)) {
    port <- CONFIG[["fronted"]][["port"]]
} else {
    port <- as.numeric(port)
}
options(shiny.port = port)
options(shiny.host = "0.0.0.0")
runApp(appDir = file.path(REICH_ROOT, "fronted"), launch.browser = F)
