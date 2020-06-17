#--------------------------------------------------------------------
# Script Name:   AnalysisDiversityVisFunctions.R
# Description:   Visualization functions for diveristy analyses
# Author:        Brandon Monier
# Created:       2020-06-17 at 11:39:52
# Last Modified: 2020-06-17 at 11:41:52
#--------------------------------------------------------------------

#--------------------------------------------------------------------
# Detailed Purpose:
#    The main purpose of this Rscript is to house functions for
#    visualizing diversity analyses of TASSEL.
#--------------------------------------------------------------------

#' @title LD visualization application
#'
#' @description Runs an interactive visualizer for an LD object
#'
#' @param ldData An LD data frame.
#'
#' @import shiny
#' @importFrom plotly layout
#' @importFrom plotly plot_ly
#' @importFrom plotly plotlyOutput
#' @importFrom plotly renderPlotly
#'
#' @export
ldApp <- function(ldData) {

    ui <- shiny::shinyUI(shiny::fluidPage(
        shiny::h4("rTASSEL - Linkage Disequilibrium"),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::numericInput(
                    inputId = "window",
                    label = "Window size",
                    value = 10
                ),
                shiny::selectInput(
                    inputId = "matVal",
                    label = "Value type",
                    choices = c(
                        "R Squared" = "R^2",
                        "P-Value" = "pDiseq",
                        "D Prime" = "DPrime"
                    )
                )
            ),
            shiny::mainPanel(
                plotly::plotlyOutput("distPlot"),
                # uiOutput("xMov"),
                # uiOutput("yMov")
                shiny::sliderInput(
                    inputId = "xMov",
                    label = "Move Window (x axis)",
                    min = 1,
                    max = 100,
                    value = 1
                ),
                shiny::sliderInput(
                    inputId = "yMov",
                    label = "Move Window (y axis)",
                    min = 1,
                    max = 100,
                    value = 1
                ),
                shiny::verbatimTextOutput("ldDebug")
            )
        )
    ))

    server <- function(input, output) {
        # output$xMov <- shiny::renderUI({
        #     shiny::sliderInput(
        #         inputId = "xMovActive",
        #         label = "Move Window (x axis)",
        #         min = 1,
        #         max = 3000,
        #         value = 1
        #     )
        # })
        # output$yMov <- shiny::renderUI({
        #     shiny::sliderInput(
        #         inputId = "yMovActive",
        #         label = "Move Window (y axis)",
        #         min = 1,
        #         max = 3000,
        #         value = 1
        #     )
        # })
        output$distPlot <- plotly::renderPlotly({
            ## Get LD matrix
            ldOut <- ldDFToShinyMat(
                ldDF = ldData,
                matVal = input$matVal,
                x_range_1 = input$xMov,
                x_range_2 = input$xMov + input$window,
                # x_range_2 = (input$xMovActive + input$window) - 1,
                y_range_1 = input$yMov,
                y_range_2 = input$yMov + input$window
                # y_range_2 = (input$yMovActive + input$window) - 1
            )

            ## Plotly metadata and parameters
            ax <- list(
                zeroline = FALSE,
                showline = FALSE,
                showticklabels = FALSE,
                showgrid = TRUE
            )

            ## The Plotly plot
            plotly::plot_ly(
                x = colnames(ldOut),
                y = rownames(ldOut),
                z = ldOut,
                type = "heatmap"
            ) %>% plotly::layout(xaxis = ax, yaxis = ax)
        })

        output$ldDebug <- shiny::renderPrint({
            debug_list <- list(
                x_range_1 = input$xMov,
                x_range_2 = input$xMov + input$window,
                y_range_1 = input$yMov,
                y_range_2 = input$yMov + input$window
            )
            cat("--- LD DEBUG ---\n")
            debug_list
        })
    }

    shiny::shinyApp(ui, server)
}


## LD dataframe to matrix converter - not exported (house keeping)

# #' @importFrom stringr str_sort
ldDFToShinyMat <- function(ldDF,
                           matVal = c("R^2", "pDiseq", "DPrime"),
                           x_range_1,
                           x_range_2,
                           y_range_1,
                           y_range_2,
                           subSet = NULL) {

    matVal <- match.arg(matVal)

    # Add new coordinates (combine chrom. and chrom. coordinate)
    ldSUB <- ldDF[, c(1:3, 7:9, 13:17)]
    ldSUB$coord1 <- paste0(ldSUB$Locus1, "_", ldSUB$Position1)
    ldSUB$coord2 <- paste0(ldSUB$Locus2, "_", ldSUB$Position2)

    # Sub matrix check
    if (!is.null(subSet)) {
        matEx <- ldSUB[1:subSet, c("coord1", "coord2", matVal)]
    } else {
        matEx <- ldSUB[, c("coord1", "coord2", matVal)]
    }

    # Subset matrix IDs
    matIDs <- unique(c(matEx$coord1, matEx$coord2))
    matIDs <- stringr::str_sort(matIDs, numeric = TRUE)

    # Create NA matrix (for population)
    mat <- matrix(data = NA, nrow = length(matIDs), ncol = length(matIDs))
    colnames(mat) <- matIDs
    rownames(mat) <- matIDs

    # Populate NA matrix with existing TASSEL calculations
    matExMat <- as.matrix(matEx)
    mat[matExMat[, 1:2]] <- as.numeric(matExMat[, 3])

    # Sub Matrix
    matSub <- mat[y_range_1:y_range_2, x_range_1:x_range_2]

    if (all(is.na(matSub))) {
        matSub[is.na(matSub)] <- 0
    }


    # Rotate and visualze
    # matCorrect <- t(apply(matSub, 2, rev))
    # matCorrect <- matSub

    # Return
    return(matSub)
    # return(dim(mat))
}


