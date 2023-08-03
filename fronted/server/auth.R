#!/usr/bin/env Rscript
library(digest)

validate_passwd <- function(passwd, hash) {
    digest(passwd, algo = "sha256", serialize = FALSE) == hash
}

do_login <- function(input, output, session) {
    auth <- reactiveValues()
    auth[["is_valid"]] <- FALSE
    auth[["username"]] <- NULL

    observeEvent(input$login, {
        user <- CONFIG[["fronted"]][["users"]][[input$username]]
        hash_passwd <- user[["hash_passwd"]]
        role <- user[["role"]]
        if (validate_passwd(input$password, hash_passwd)) {
            auth[["is_valid"]] <- TRUE
            auth[["username"]] <- input$username
            auth[["role"]] <- role
            showNotification("Login successful!", type = "message")
        } else {
            showModal(modalDialog(
                title = "Invalid Credentials",
                "The username or password you entered is incorrect.",
                easyClose = TRUE,
                footer = NULL
            ))
        }
    })

    return(auth)
}
