
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Import library                          ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(shiny)
library(plotly)
library(bslib)
library(knitr)
library(markdown)
library(BiocParallel)
library(ggplot2)
library(GCIMS)
library(shinyFiles)
library(shinycssloaders)
library(prompter)
library(yaml)
library(fs)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Import function                         ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' Image Compression
#'
#' Reduce the precision of an image by making n*n mean through the image.
#' Each n*n squares of pixels are replaced by 1 mean pixel. It divides by
#' n the number of pixels in the image.
#'
#' @param M matrix of the intensity. A numeric or logical matrix containing
#' the values to be plotted cf. image().
#' @param x abscissa, locations of grid lines at which the values in z are
#' measured cf. image().
#' @param y ordinate, locations of grid lines at which the values in z are
#' measured cf. image().
#' @param n division of the precision.
#'
#' @return M_comp which is the matrix M with its precision
#' reduced by n. x_comp and y_comp, linear vectors containing the
#' corresponding compressed abscissa and ordinary.
#' @export
#'
#' @examples
#' image_comp(
#'   matrix(1:36, nrow = 6, ncol = 6, byrow = TRUE),
#'   x = 1:6,
#'   y = 1:6,
#'   n = 2
#' )
image_comp<- function(M,x,y,n){

  size<-dim(M)# M matrix dimensions

  ##Matrix compression
  M_comp<-sapply(seq(1,size[2]-(n-1), by=n),function(j) #while we cross M lines
  {
    sapply(seq(1,size[1]-(n-1), by=n),function(i) #cross M columns
    {
      m<-mean(M[i:(i+n-1),j:(j+n-1)]) #do the div*div mean
    })
  })#sapply returns the matrix directly

  ##Abscissa compression
  x_comp<-sapply(seq(1,size[1]-(n-1), by=n),function(k) #while we cross M lines
  {
    m<-mean(x[k:(k+n-1)])
  })
  ##Ordinate compression
  y_comp<-sapply(seq(1,size[2]-(n-1), by=n),function(k) #while we cross M lines
  {
    m<-mean(y[k:(k+n-1)])
  })

  return(list("M_comp"=M_comp,"x_comp"=x_comp,"y_comp"=y_comp))
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Define UI                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fluidPage(

  theme = bs_theme(fg = "rgb(58, 162, 157)", font_scale = NULL,
                   bootswatch = "materia", bg = "rgb(255, 255, 255)",
                   version=5
  ),#bootstrap theme
  #App theme  ----

  #Button animation
  use_prompt(),

  #Nav Bar  ----
  navbarPage("GC-IMS GUI",

             ##~~~~~~~~~~~~~~~~~~~~~~~~
             ##  ~ First tab : Data Preparation ----
             ##~~~~~~~~~~~~~~~~~~~~~~~~
             tabPanel("Data Preparation",

                      # Sidebar layout with input and output definitions--------
                      sidebarLayout(

                        # Sidebar panel for inputs
                        sidebarPanel(
                          h1("Data Preparation"),

                          br(),

                          p("First, you need to load an annotations file which is
                          based on your samples. It contains the id, the name
                          and the metadata of your  samples. If you don't have
                          any annotations file, create one."),

                          h2("Create an annotations file"),

                          fluidRow(
                            column(7,
                                   verbatimTextOutput("samplesdir"),
                            ),
                            column(2,
                                   shinyDirButton("samples_dir",
                                                  "Folder",
                                                  "Please select a folder"),
                            ),
                          ),#fluidrow

                          br(),

                          shinySaveButton("annotations_dir",
                                          "Save annotations file",
                                          "Save file as...",
                                          filetype = list(excel="xlsx"),
                                          viewtype = "icon"),

                          fluidRow(
                            add_prompt(
                              fileInput("browse_annotations",
                                        label = h2("Choose your annotations file")
                              ),
                              position = "bottom",
                              message = ".csv or .xlsx"
                            ),#add_prompt
                          ),#fluidRow

                        ),#sidebar panel

                        # Main panel for displaying outputs---------------------
                        mainPanel(

                          tabsetPanel(

                            tabPanel("Conversion",
                                     br(),
                                     h6("In this window you can convert your raw
                                     data into .rds files."),
                                     br(),

                                     #Select the folder containing data to be converted
                                     fluidRow(
                                       column(8,
                                              verbatimTextOutput("rawdir"),
                                       ),#column
                                       column(4,
                                              shinyDirButton("raw_dir",
                                                             "Raw folder",
                                                             "Please select a folder"),

                                       ),#column
                                     ),#fluidRow

                                     #Select a folder to store converted data
                                     fluidRow(
                                       column(8,
                                              verbatimTextOutput("convertdir"),
                                       ),#column
                                       column(4,
                                              shinyDirButton("convert_dir",
                                                             "Save in",
                                                             "Please select a folder"),
                                       ),#column
                                     ),#fluidRow

                            ),#tabpanel conversion
                          )#tabsetpanel
                        )#mainPanel
                      )#sidebarLayout
             ),#tabPanel

             ##~~~~~~~~~~~~~~~~~~~~~~~~
             ##  ~ Second tab : Plot ----
             ##~~~~~~~~~~~~~~~~~~~~~~~~
             tabPanel("Plot",

                      # Sidebar layout with input and output definitions ----
                      sidebarLayout(

                        # Sidebar panel for inputs ----
                        sidebarPanel(
                          h1("Dynamic plot"),
                          p("Visualise and interact with the data."),

                          #Select a sample to be plotted
                          selectInput('filestoplot_name', h5('Choose the file you want to plot'),
                                      choices=NULL, selectize=TRUE),

                          h2("Parameters"),

                          fluidRow(
                            column(12,
                                   fileInput("browse_parameters",
                                             label = h5("Load your parameters")
                                   ),#fileInput
                            ),#column
                          ),#fluidRow

                          tabsetPanel(
                            tabPanel("Reshape & Precision",
                                     br(),
                                     fluidRow(
                                       column(6,
                                              sliderInput("dt_reshape", "Drift time range :",
                                                          min = 0, max = 30,
                                                          value = c(1,30)),
                                       ),#column
                                       column(6,
                                              sliderInput("rt_reshape", "Retention time range :",
                                                          min = 0, max = 2000,
                                                          value = c(0,2000)),
                                       ),#column
                                     ),#fluidRow
                                     sliderInput("precision", p("Precision :"),
                                                 min = 0, max = 10,
                                                 value = 4),
                            ),#tabPanel

                            tabPanel("Savitzky-Golay filter",

                                     br(),

                                     p("Drift time"),

                                     fluidRow(
                                       column(2,
                                              radioButtons("SG_dt_order", label = "order :",
                                                           choices = list( 1, 2, 3),
                                                           selected = 2,
                                                           inline=T),
                                       ),#column
                                       column(5,
                                              sliderInput("SG_dt_size", "length (ms) :",
                                                          min = 0, max = 1,
                                                          value = 0.14),
                                       ),#column
                                       column(2,
                                              br(),
                                              actionButton("SG_dt_apply", "apply")),#column
                                     ),#fluidRow

                                     p("Retention time"),

                                     fluidRow(
                                       column(2,
                                              radioButtons("SG_rt_order", label = "order :",
                                                           choices = list( 1, 2, 3),
                                                           selected = 2,
                                                           inline=T),
                                       ),#column
                                       column(5,
                                              sliderInput("SG_rt_size", "length (s) :",
                                                          min = 3, max = 30,
                                                          value = 3),
                                       ),#column
                                       column(2,
                                              br(),
                                              actionButton("SG_rt_apply", "apply")),#column
                                     ),#fluidRow

                                     actionButton("SG_apply", "apply to the whole sample"),
                            ),#tabPanel


                            tabPanel("File",
                                     br(),
                                     p("Save your parameters"),

                                     shinySaveButton("yaml_dir",
                                                     "Save",
                                                     "Save file as...",
                                                     filetype = list(excel="yaml"),
                                                     viewtype = "icon"),
                            )#tabPanel
                          ),#tabsetPanel

                          width=3
                        ),#sidebarPanel

                        # Main panel for displaying outputs ----
                        mainPanel(
                          fluidRow(column(6,withSpinner(plotlyOutput("plot_RT"))),
                                   column(6,withSpinner(plotlyOutput("mainplot")))),
                          fluidRow(column(6, verbatimTextOutput("rt_dt_click")),
                                   column(6,withSpinner(plotlyOutput("plot_DT"))))
                        )#mainPanel
                      )#sidebarLayout
             ),#tabPanel

             ##~~~~~~~~~~~~~~~~~~~~~~~~
             ##  ~ Third tab : Pre Processing ----
             ##~~~~~~~~~~~~~~~~~~~~~~~~
             #tabPanel("Pre Processing",paste0("Pre Processing")),

             ##~~~~~~~~~~~~~~~~~~~~~~~~
             ##  ~ Fourth tab : Help ----
             ##~~~~~~~~~~~~~~~~~~~~~~~~
             tabPanel("Help",
                      paste0("help"),
                      #img(src = "logo.png", style = "float: left; width: 120px; margin-right: 10px; margin-top: 5px"),
             )
  )
)

