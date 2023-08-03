#!/usr/bin/env Rscript


server <- function(input, output, session) {
    source("server/auth.R", local = TRUE)

    auth <- do_login(input, output, session)
    output$main_ui <- renderUI({
        if (auth[["is_valid"]]) {
            source("ui/main.R")
            main_ui
        } else {
            source("ui/login.R")
            login_ui
        }
    })
}
