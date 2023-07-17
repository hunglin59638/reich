#!/usr/bin/env Rscript
library(digest)

validate_passwd <- function(passwd, hash) {
    digest(passwd, algo = "sha256", serialize = FALSE) == hash
}
