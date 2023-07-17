#!/usr/bin/env Rscript

source("server/auth.R")
hash_passwd <- "8c6976e5b5410415bde908bd4dee15dfb167a9c873fc4bb8a81f6f2ab448a918"
server <- function(input, output, session) {
    user <- reactiveValues()
    user[["is_valid"]] <- FALSE
    user[["username"]] <- NULL

    observeEvent(input$login, {
        if (validate_passwd(input$password, hash_passwd)) {
            user[["is_valid"]] <- TRUE
            user[["username"]] <- input$username
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

    output$main_ui <- renderUI({
        if (user[["is_valid"]]) {
            source("ui/main.R")
            main_ui
        } else {
            source("ui/login.R")
            login_ui
        }
    })
}
