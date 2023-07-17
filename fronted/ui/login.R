#!/usr/bin/env Rscript

login_ui <- fluidPage(
    tags$style(
        HTML(
            "
            .center_box {
                background-color: lightgray;
                border: 1px solid gray;
                border-radius: 20px;
                padding: 20px;
                text-align: center;
            }
            "
        )
    ),
    absolutePanel(
        top = "30%",
        left = "40%",
        style = "
        background-color: lightgray;
        border-radius: 20px;
        ",
        height = 250,
        width = 400,
        wellPanel(
            style = "background-color: lightgray;",
            h3("Reich", class = "center_box"),
            br(),
            textInput("username", "Username: "),
            br(),
            passwordInput("password", "Password: "),
            br(),
            actionButton("login", "Login")
        )
    )
)
