
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
##                                  Define server                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

function(input, output, session) {

  ##~~~~~~~~~~~~~~~~~~~~~~~~
  ##  ~ Data preparation ----
  ##~~~~~~~~~~~~~~~~~~~~~~~~

  #Home path
  volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), getVolumes()())

  #Create an annotations file---------------------------------------------------

  #Get the directory of the folder containing the samples you'll work with

  shinyDirChoose(input,
                 "samples_dir",
                 roots = volumes,
                 session = session,
                 restrictions = system.file(package = "base"),
                 allowDirCreate = TRUE)

  output$samplesdir <- renderText({#display the directory chosen
    if (is.integer(input$samples_dir)) {
      paste0("Select a folder with your raw samples")
    } else {
      parseDirPath(volumes, input$samples_dir)
    }
  })

  samples_directory<-reactive({value=parseDirPath(volumes, input$samples_dir)})

  #Annotations directory

  shinyFileSave(input, "annotations_dir",
                roots = volumes,
                session = session,
                restrictions = system.file(package = "base"))

  annotations_info<-reactive({value=parseSavePath(volumes, input$annotations_dir)})

  annotations_directory<-reactive({value=annotations_info()$datapath})

  observeEvent(annotations_directory(), {
    if (length(annotations_directory())==!0){
      # Your samples match this extension:
      samples_ext <- "*.mea.gz"
      filenames <- list.files(samples_directory(), pattern = utils::glob2rx(samples_ext),
                              recursive = TRUE, include.dirs = FALSE)
      #
      annotations_template <- data.frame(
        SampleID = tools::file_path_sans_ext(filenames, compression = TRUE),
        FileName = filenames
      )
      # You may need to run: install.packages("writexl")
      writexl::write_xlsx(annotations_template, annotations_directory())
    }
  })


  # Disable parallellization: (Useful for better error reporting)--------------

  register(SerialParam(progressbar = TRUE), default = TRUE) #QUESTION KESAKO

  #Browse your annotations file-------------------------------------------------

  annotationsfile<-reactive({value=input$browse_annotations})

  #make the annotations variable global
  annotations <- NULL
  makeReactiveBinding("annotations")

  observeEvent(annotationsfile(), {
    #save the content of your annotations file in a df variable named "annotations"
    if (tolower(tools::file_ext(annotationsfile()$name))=="xlsx"){
      print("the extention of the file is xlsx")
      annotations <<- readxl::read_excel(annotationsfile()$datapath)
    }

    if(tolower(tools::file_ext(annotationsfile()$name))=="csv"){
      print("the extention of the file is csv")
      annotations <<- readr::read_csv(annotationsfile()$datapath, show_col_types = FALSE)
    }
    print(annotations)
  })


  #Conversion-------------------------------------------------------------------

  #Folder containing raw data to be converted
  shinyDirChoose(input,
                 "raw_dir",
                 roots = volumes,
                 session = session,
                 restrictions = system.file(package = "base"),
                 allowDirCreate = TRUE)

  output$rawdir <- renderText({
    if (is.integer(input$raw_dir)) {
      paste0("Select a folder with raw samples to be loaded")
    } else {
      parseDirPath(volumes, input$raw_dir)
    }
  })

  raw_directory<-reactive({value=parseDirPath(volumes, input$raw_dir)})

  #Folder to save converted data
  shinyDirChoose(input,
                 "convert_dir",
                 roots = volumes,
                 session = session,
                 restrictions = system.file(package = "base"),
                 allowDirCreate = TRUE)

  output$convertdir<- renderText({
    if (is.integer(input$convert_dir)) {
      paste0("Select a folder to store converted data")
    } else {
      parseDirPath(volumes, input$convert_dir)
    }
  })

  convert_directory<-reactive({value=parseDirPath(volumes, input$convert_dir)})

  #make the dataset variable global
  dataset <- NULL
  makeReactiveBinding("dataset")

  observeEvent(convert_directory(), {
    if (length(convert_directory())==!0){
      #create GCIMSdataset object
      dataset<<-GCIMSDataset(annotations,
                             base_dir = raw_directory(),
                             scratch_dir = convert_directory())
      getGCIMSSample(dataset, sample = annotations$SampleID[1])
      print(dataset)
    }
  })#observeEvent


  ##~~~~~~~~~~~~~~~~~~~~~~~~
  ##  ~ Plot ----
  ##~~~~~~~~~~~~~~~~~~~~~~~~

  #Choose a file you want to plot ----------------------------------------------

  #Update the selection bar
  observeEvent(annotationsfile(), {
    choices <- annotations$SampleID
    updateSelectInput(inputId = "filestoplot_name", choices = choices)
  })

  #Select a file to plot
  filestoplot_select<-reactive({value=input$filestoplot_name})

  #Load data to plot -----------------------------------------------------------

  untreated_sample<-reactive({value=getGCIMSSample(dataset, sample = filestoplot_select())})
  untreated_sample_intensity<-reactive({value=untreated_sample()@data})
  untreated_sample_dt<-reactive({value=untreated_sample()@drift_time})
  untreated_sample_rt<-reactive({value=untreated_sample()@retention_time})

  #Save parameters (optional)----------------------------------------------

  parameters<-reactive({value=data.frame(precision=input$precision,
                                         SG_dt_order=input$SG_dt_order,
                                         SG_dt_size=input$SG_dt_size,
                                         SG_rt_order=input$SG_rt_order,
                                         SG_rt_size=input$SG_rt_size,
                                         dt_range_min=input$dt_reshape[1],
                                         dt_range_max=input$dt_reshape[2],
                                         rt_range_min=input$rt_reshape[1],
                                         rt_range_max=input$rt_reshape[2])})



  shinyFileSave(input, "yaml_dir",
                roots = volumes,
                session = session,
                restrictions = system.file(package = "base"))

  yaml_info<-reactive({value=parseSavePath(volumes,input$yaml_dir)})
  yaml_directory<-reactive({value=yaml_info()$datapath})

  observeEvent(yaml_directory(),{
    if (length(yaml_directory())==!0){
      print("hey!")
      write_yaml(parameters(), yaml_directory())
    }
  })

  #Load your parameters (optional)----------------------------------------------

  parametersfile<-reactive({value=input$browse_parameters})

  observeEvent(input$browse_parameters,{
    parameters<-read_yaml(file =parametersfile()$datapath)
    print(parameters)
    updateSliderInput(inputId = "dt_reshape", value=c(parameters$dt_range_min,parameters$dt_range_max))
    updateSliderInput(inputId = "rt_reshape", value=c(parameters$rt_range_min,parameters$rt_range_max))
    updateSliderInput(inputId = "precision", value=parameters$precision)
    updateSliderInput(inputId = "SG_dt_size", value=parameters$SG_dt_size)
    updateSliderInput(inputId = "SG_rt_size", value=parameters$SG_rt_size)
    updateRadioButtons(inputId = "SG_dt_order",selected=parameters$SG_dt_order)
    updateRadioButtons(inputId = "SG_rt_order",selected=parameters$SG_rt_order)
  })


  #Reshape----------------------------------------------------------------------

  #Get range
  range_dt<-reactive({value=input$dt_reshape})
  range_rt<-reactive({value=input$rt_reshape})

  #Initialization and make variables global
  sampleafterfilter<-reactive({value=getGCIMSSample(dataset, sample = filestoplot_select())})
  sampleafterfilter_intensity<-reactive({value=sampleafterfilter()@data})
  sampleafterfilter_dt<-reactive({value=sampleafterfilter()@drift_time})
  sampleafterfilter_rt<-reactive({value=sampleafterfilter()@retention_time})

  makeReactiveBinding("sampleafterfilter")
  makeReactiveBinding("sampleafterfilter_intensity")
  makeReactiveBinding("sampleafterfilter_dt")
  makeReactiveBinding("sampleafterfilter_rt")

  #Update when the slider moves
  observeEvent(range_dt(),{
    observeEvent(dataset,{#to avoid crashing when dataset doesn't exist yet
      sampleafterfilter<<-reactive({value=filterDt(untreated_sample(), dt = range_dt())}) # in ms
      sampleafterfilter_intensity<<-reactive({value=sampleafterfilter()@data})
      sampleafterfilter_dt<<-reactive({value=sampleafterfilter()@drift_time})
      sampleafterfilter_rt<<-reactive({value=sampleafterfilter()@retention_time})
    })
  })

  observeEvent(range_rt(),{
    observeEvent(dataset,{
      sampleafterfilter<<-reactive({value=filterRt(untreated_sample(), rt = range_rt())}) # in s
      sampleafterfilter_intensity<<-reactive({value=sampleafterfilter()@data})
      sampleafterfilter_dt<<-reactive({value=sampleafterfilter()@drift_time})
      sampleafterfilter_rt<<-reactive({value=sampleafterfilter()@retention_time})
    })
  })

  #SG filter--------------------------------------------------------------------

  #Drift Time
  SG_dt_order<-reactive({value=input$SG_dt_order})
  SG_dt_size<-reactive({value=input$SG_dt_size})

  #Get one intensity line
  one_ims_spec <- reactive({value=getIMS(sampleafterfilter(), rt_range = clicked_rt())})
  one_ims_smoothed <- reactive({value=smooth(one_ims_spec(), dt_length_ms = SG_dt_size(), dt_order = strtoi(SG_dt_order()))})

  #Make button variable global
  dt_smooth_apply<-FALSE
  makeReactiveBinding("dt_smooth_apply")

  observeEvent(input$SG_dt_apply,{
    dt_smooth_apply<<-!dt_smooth_apply
  })

  #Retention Time
  SG_rt_order<-reactive({value=input$SG_rt_order})
  SG_rt_size<-reactive({value=input$SG_rt_size})

  #Get one intensity column
  one_chrom <- reactive({value=getEIC(sampleafterfilter(), dt_range = clicked_dt())})
  one_chrom_smoothed <- reactive({value=smooth(one_chrom(), rt_length_s = SG_rt_size(), rt_order = strtoi(SG_rt_order()))})

  #Make button variable global
  rt_smooth_apply<-FALSE
  makeReactiveBinding("rt_smooth_apply")

  observeEvent(input$SG_rt_apply,{
    rt_smooth_apply <<- !rt_smooth_apply
  })

  #Precision and display--------------------------------------------------------

  #Get precision slider value
  precisionvalue<- reactive({value=input$precision})

  #Function that reduces the precision of the image
  sample<-reactive({value=image_comp(sampleafterfilter_intensity(),
                                     sampleafterfilter_dt(),
                                     sampleafterfilter_rt(),
                                     11-precisionvalue())
  })
  sample_dt<-reactive({value=sample()$x_comp})
  sample_rt<-reactive({value=sample()$y_comp})
  sample_intensity<-reactive({value=sample()$M_comp})

  #sharpen the intensity to have a better display
  #(I transposed because it was different formalism)
  sample_display<-reactive({value=t(sign(sample_intensity())* abs(sample_intensity())^(1/3))})

  #Plots section----------------------------------------------------------------

  output$rt_dt_click<-renderPrint({
    paste0("dt=",clicked_dt(),"ms rt=",clicked_rt(),"s")
  })

  #Main Plot
  output$mainplot <- renderPlotly({

    if (is.null(dtmin())){#if there is no zoom
      dtmin<-min(sample_dt())
      dtmax<-max(sample_dt())
      rtmin<-min(sample_rt())
      rtmax<-max(sample_rt())
    }
    else{#if there is a zoom
      dtmin<-dtmin()
      dtmax<-dtmax()
      rtmin<-rtmin()
      rtmax<-rtmax()
    }

    fig <- plot_ly(
      x = sample_dt(),
      y = sample_rt(),
      z = sample_display(),
      type = "heatmap",
      showscale = FALSE,
      source = "source")

    fig=layout( fig,
                xaxis = list(
                  range = list(dtmin,dtmax)),
                yaxis = list(
                  range = list(rtmin,rtmax)))
  })

  #Zoom info
  zoom<-reactive({value=event_data("plotly_relayout",source = "source")})
  dtmin<-reactive({value=zoom()$`xaxis.range[0]`})
  dtmax<-reactive({value=zoom()$`xaxis.range[1]`})
  rtmin<-reactive({value=zoom()$`yaxis.range[0]`})
  rtmax<-reactive({value=zoom()$`yaxis.range[1]`})
  index_rtmin<-reactive({value=which.min(abs(sampleafterfilter_rt()-rtmin()))})
  index_rtmax<-reactive({value=which.min(abs(sampleafterfilter_rt()-rtmax()))})
  index_dtmin<-reactive({value=which.min(abs(sampleafterfilter_dt()-dtmin()))})
  index_dtmax<-reactive({value=which.min(abs(sampleafterfilter_dt()-dtmax()))})

  #Click info
  clicked <- reactive({value=event_data("plotly_click",source = "source")})
  clicked_rt<-reactive({value=clicked()$y})
  clicked_dt<-reactive({value=clicked()$x})

  #RT Plot
  output$plot_RT <- renderPlotly({

    #Instance
    dt_ms <- sampleafterfilter_dt()
    rt_s <- sampleafterfilter_rt()
    intens <- sampleafterfilter_intensity()
    index_clicked_dt<-which.min(abs(dt_ms-clicked_dt()))

    #Default opening value
    if(is.null(clicked_dt())) index_clicked_dt<-1

    #Display

    if(rt_smooth_apply==TRUE) {

      chrom <- one_chrom()
      if (is.null(rtmin())){
        fig2 <-plot_ly(
          x=intensity(chrom),
          y=rtime(chrom),
          type = 'scatter',
          mode = 'lines'
        )

        fig2<-add_trace(
          fig2,
          type = 'scatter',
          mode = 'line',
          x=one_chrom_smoothed()@intensity,
          y=one_chrom_smoothed()@retention_time
        )

      } else{
        fig2 <-plot_ly(
          x=one_chrom()@intensity[index_rtmin():index_rtmax()],
          y=one_chrom()@retention_time[index_rtmin():index_rtmax()],
          type = 'scatter',
          mode = 'lines'
        )
        fig2<-add_trace(
          fig2,
          type = 'scatter',
          mode = 'line',
          x=one_chrom_smoothed()@intensity[index_rtmin():index_rtmax()],
          y=one_chrom_smoothed()@retention_time[index_rtmin():index_rtmax()]
        )

      }

    }else{

      if (is.null(rtmin())){
        fig2 <-plot_ly(
          x=intens[index_clicked_dt,],
          y=rt_s,
          type = 'scatter',
          mode = 'lines'
        )
      } else{
        fig2 <-plot_ly(
          x=intens[index_clicked_dt,index_rtmin():index_rtmax()],
          y=rt_s[index_rtmin():index_rtmax()],
          type = 'scatter',
          mode = 'lines'
        )

      }

    }

    fig2=layout( fig2,
                 xaxis = list(
                   title = 'Intensity'),
                 yaxis = list(
                   title = 'Retention Time (s)'),
                 showlegend = FALSE
    )
  })

  #DT Plot
  output$plot_DT <- renderPlotly({

    #Instance
    dt_ms <- sampleafterfilter_dt()
    rt_s <- sampleafterfilter_rt()
    intens <- sampleafterfilter_intensity()
    index_clicked_rt<-which.min(abs(rt_s-clicked_rt()))

    #Default opening value
    if(is.null(clicked_rt())) index_clicked_rt<-1

    #Display

    if(dt_smooth_apply==TRUE){

      if (is.null(dtmin())){
        fig3 <-plot_ly(
          x=one_ims_spec()@drift_time,
          y=one_ims_spec()@intensity,
          type = 'scatter',
          mode = 'lines')

        fig3<-add_trace(
          fig3,
          type = 'scatter',
          mode = 'line',
          x=one_ims_smoothed()@drift_time,
          y=one_ims_smoothed()@intensity
        )

      } else{
        fig3 <-plot_ly(
          x=one_ims_spec()@drift_time[index_dtmin():index_dtmax()],
          y=one_ims_spec()@intensity[index_dtmin():index_dtmax()],
          type = 'scatter',
          mode = 'lines')

        fig3<-add_trace(
          fig3,
          type = 'scatter',
          mode = 'line',
          x=one_ims_smoothed()@drift_time[index_dtmin():index_dtmax()],
          y=one_ims_smoothed()@intensity[index_dtmin():index_dtmax()]
        )
      }

    }else{

      if (is.null(dtmin())){
        fig3 <-plot_ly(
          x=dt_ms,
          y=intens[,index_clicked_rt],
          type = 'scatter',
          mode = 'lines')

      } else{
        fig3 <-plot_ly(
          x=dt_ms[index_dtmin():index_dtmax()],
          y=intens[index_dtmin():index_dtmax(),index_clicked_rt],
          type = 'scatter',
          mode = 'lines')
      }
    }

    fig3=layout( fig3,
                 xaxis = list(
                   title = 'Drift Time (ms)'),
                 yaxis = list(
                   title = 'Intensity'),
                 showlegend = FALSE
    )

  })

}


