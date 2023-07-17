#!/usr/bin/env Rscript

library(shiny)
library(glue)

data_panel <- tabPanel(
    "My Data",
    br(),
)

upload_panel <- tabPanel(
    "Upload",
    fluidPage(
        fluidRow(
            column(2),
            column(
                7,
                h3("Upload Samples"),
                hr(),
                selectInput("project", "Select Project", choice = c(), width = 1200),
                actionLink("add_project", "+ CREATE PROJECT"),
                br(),
                h3("Analysis Type"),
                radioButtons("analysis_type", br(),
                    choiceNames = c("Metagenomics", "SARS-CoV-2 Consensus Genome"),
                    choiceValues = c("metagenomics", "covid")
                ),
                h3("Upload Files"),
                tabsetPanel(
                    tabPanel(
                        "Upload from Your Computer",
                        fileInput("fq_files", br(),
                            accept = c(".fq", ".fastq", ".fq.gz", ".fastq.gz"),
                            multiple = TRUE,
                            width = 1200,
                            placeholder = "Drag and drop your files here"
                        ),
                    ),
                    tabPanel("Upload from Basespace")
                ),
                fluidRow(
                    column(
                        1,
                        actionButton("continue", "Continue")
                    ),
                    column(
                        1,
                        actionButton("concel", "Cancel")
                    ),
                    column(10),
                )
            ),
            column(3)
        )
    )
)
bottom_line_cls <- "bottom_line"

header <- conditionalPanel(
    condition = "input.tabs == 'My Data'",
    fluidPage(
        tags$style(glue(".{bottom_line_cls} {{ border-bottom: 1px solid lightgrey; }}")),
        fluidRow(
            column(
                2,
                fluidRow(
                    column(2, icon("home")),
                    column(10, textInput("search", NULL, placeholder = "Search My Data..."))
                )
            ),
            column(
                2,
                fluidRow(
                    column(3, actionLink("activate_projects", "Projects")),
                    column(3, actionLink("activate_samples", "Samples")),
                    column(3, actionLink("activate_visualizations", "Visualizations"))
                )
            ),
        ),
        fluidRow(
            div(class = bottom_line_cls, column(12))
        )
    )
)

main_ui <- navbarPage(
    id = "tabs",
    title = "Reich",
    position = "static-top",
    inverse = FALSE,
    selected = "My Data",
    header = header,
    data_panel,
    upload_panel,
)
