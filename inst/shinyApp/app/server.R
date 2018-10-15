
# library("shinyFiles")
library("shiny")
library("shinyjs")

## install all packages and check
list.of.packages <- c("ggplot2", "Rcpp","cluster","igraph","plyr","reshape","scales",
  "grDevices","parallel","jsonlite","doParallel","shiny","shinydashboard","shinyjs","gtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
list.of.biocondu <- c("flowCore","ggcyto")
new.packages <- list.of.biocondu[!(list.of.biocondu %in% installed.packages()[,"Package"])]
if(length(new.packages)){source("https://bioconductor.org/biocLite.R"); biocLite("flowCore"); biocLite("ggcyto");}

options(expressions = 5e5,shiny.maxRequestSize = 20 * 1024 ^ 3) #increase max file upload size to 20 Gb

shinyServer(function(input, output, session) {
  session$onSessionEnded(stopApp)
  
  progress <- Progress$new()
  progress$set(message="Load library", value=0.33)
  library("Rcpp")
  library("cluster")
  library("flowCore")
  library("ggplot2")
  library("igraph")
  library("plyr")
  library("reshape")
  progress$set(message="Load library", value=0.66)
  library("scales")
  library("grDevices")
  library("parallel")
  library("jsonlite")
  library("doParallel")
  progress$set(message="Read Source", value=1)
  source("functions.R")
  sourceCpp("forceatlas2.cpp")
  progress$close()
  
  listObject <- reactiveValues(
    ##GLOBAL
    # outputDirectory = paste0(getwd(), "/Data/TEST 3"),
    # outputDirectory = "C:/Users/Cmatteoli/Documents/SCAFFOLD/Data/test/generated vs IMPC/maps",
    outputDirectory = getwd(),
    inputDirectory = getwd(),
    clustering.groups = NULL,
    files.id = NULL,
    origin.parameters = NULL,
    params.ori = NULL,
    param.ref = NULL,

    ##CLUSTERING
    doneComp = NULL,
    flow.frames = NULL,
    over.clustering = NULL,
    flow.frames.enrich = NULL,
    datapaths = NULL,
    transform.data = NULL,
    cofactor = NULL,
    clustering_start = "NO",
    clustering_check = "NO",

    ##ANALYSIS
    scaffoldCoupleFiles = NULL,
    analysisFiles = list(),
    couplesid = NULL,
    inputFiles = NULL,

    a.cofactor = NULL,
    a.transform.data = NULL,

    loaded.rdata = list(),
    clustered.tables = NULL,
    clustered.txt = NULL,
    gated.flow.frames = NULL,
    gated.datapaths = NULL,
    gated.files.id = NULL,

    analysis_start = "NO",
    analysis_check = "NO",

    ##MAPPING
    file.scaffold = NULL,
    scaffold.data = NULL,

    ##MAP DATASET
    d.file.scaffold.dataset = NULL,
    d.file.clustered.ref = NULL,
    d.files.clustered.input = NULL,
    d.files.rdata = list(),
    d.files.clustered.tables = NULL,
    d.clustered.txt = NULL,
    d.files.clustered.dataset = list(),
    d.files.clustered.id = NULL,

    mapping_start = "NO",
    mapping_check = "NO"
  )

########################################################################################
############################ MAKE GROUPS FILES CLUS ####################################

  output$boxInput <- renderUI({
    if(input$tab != "mapping"){
      sep <- list("\t",",",".",";"," ")
      names(sep) <- c("tab",",",".",";","space")
      box(title="Choose FCS or TxT File(s)",collapsible = TRUE,status = "primary", solidHeader = TRUE, width = 250,
        fileInput("fileInput", "",
                multiple = TRUE,
                accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".fcs",".csv"
                        )
        ),
        fluidRow(
          column(6, 
            if(input$tab == "multiGroups"){
              selectInput("separator","CSV Separator",choices = sep, selected=1)
            }),
          column(6, actionButton("refreshInput","Refresh Input"))
        )     
      )
    }
  })

  observeEvent(input$refreshInput,{
    listObject$flow.frames <- NULL
    listObject$flow.frames.enrich <- NULL
    listObject$over.clustering <- NULL
  })

  output$boxInputTable <- renderUI({
    if(is.null(listObject$flow.frames) && is.null(listObject$over.clustering)) {return(NULL)}
    if(input$tab != "mapping"){
      if(is.null(listObject$over.clustering)){
        files <- names(listObject$flow.frames)
      } else if(is.null(listObject$flow.frames)){
        files <- names(listObject$over.clustering)
      }else{
        files <- c("ERROR NO SAME FORMAT FILES !!",names(listObject$over.clustering),names(listObject$flow.frames))
      }
      temp.files <- as.matrix(files)
      colnames(temp.files) <- "Names"
      box(title="List Files",collapsible = TRUE,status = "primary", solidHeader = TRUE, width = 250,
        renderTable({temp.files})
      )
    }
  })

  observeEvent(input$fileInput,{
    progress <- Progress$new()
    progress$set(message = "Read Upload file", value = 0)
    filesname <-c(names(listObject$flow.frames),names(listObject$over.clustering),as.vector(input$fileInput$name))

    if(length(grep(".fcs$",filesname))== length(filesname)) {
      progress$set(message="Read FCS", value=0)
      listObject$datapaths <- input$fileInput$datapath

      i <- 0
      new.flow.frames <-lapply(as.vector(listObject$datapaths), function(x) {
        i <<- i+1
        progress$set(message = paste0("Reading file ", i, "/", length(listObject$datapaths), "."), value=i/length(as.vector(listObject$datapaths)))
        return(read.FCS(x, emptyValue = FALSE))
      })

      listObject$flow.frames  <- c(listObject$flow.frames, new.flow.frames)
      names(listObject$flow.frames) <- filesname
      output$warning <- renderUI({
        textOutput({return("")})
      })

      labels <- pData(parameters(listObject$flow.frames[[1]]))[,2]
      params.names <- pData(parameters(listObject$flow.frames[[1]]))[,1]
      labels[which(is.na(labels))] <- colnames(listObject$flow.frames[[1]])[c(which(is.na(labels)))]
      labels[which(labels=="<NA>")] <- colnames(listObject$flow.frames[[1]])[c(which(labels=="<NA>"))]
      names(params.names) <- labels
      listObject$params.ori <- params.names

    }else if(length(grep(".txt$",filesname))==length(filesname) || length(grep(".csv$",filesname))==length(filesname)) {
      progress$set(message="Read FCS", value=0)
      listObject$datapaths <- input$fileInput$datapath

      i <- 0
      new.over.clustering <-lapply(as.vector(listObject$datapaths), function(x) {
        i <<- i+1
        progress$set(value=i/length(as.vector(listObject$datapaths)))
        return(read.csv(x, sep=input$separator, check.names = F))
      })

      listObject$over.clustering  <- c(listObject$over.clustering, new.over.clustering)
      names(listObject$over.clustering) <- filesname
      output$warning <- renderUI({
        textOutput({return("")})
      })
    } else {
      output$warning <- renderUI({
        textOutput({return("Error input format !")})
      })
    }
    progress$close()
  })

  observeEvent(input$separator,{
    if(is.null(input$fileInput)) return(NULL)
    progress <- Progress$new()
    progress$set(message = "Read Upload file", value = 0)
    filesname <-names(listObject$over.clustering)
    i <- 0
    new.over.clustering <-lapply(as.vector(listObject$datapaths), function(x) {
      i <<- i+1
      progress$set(value=i/length(as.vector(listObject$datapaths)))
      return(read.csv(x, sep=input$separator, check.names = F))
    })
    listObject$over.clustering  <- new.over.clustering
    names(listObject$over.clustering) <- filesname
    progress$close()
  })

  output$boxOutput <- renderUI({
    if(is.null(listObject$flow.frames) && is.null(listObject$over.clustering)) return(NULL)
    if(input$tab != "mapping"){
      output$output_dir <- renderText(return(listObject$outputDirectory))
      box(title="Select Output",collapsible = TRUE,status = "primary", solidHeader = TRUE, width = 250,
        actionButton("saveClustering", "Saving Folder"), br(), br(),
        verbatimTextOutput("output_dir")
      )
    }
  })

  observeEvent(input$saveClustering,{
    listObject$outputDirectory <- chooseDir()
    if(is.na(listObject$outputDirectory)){
      listObject$outputDirectory <- NULL
    }
  })  

########################################################################################
############################ MAKE GROUPS FILES CLUS ####################################
  ## DELETE BOX 

  output$flow_frames <- renderTable({
    if(!is.null(listObject$flow.frames)) {
      table1 <- as.matrix(names(listObject$flow.frames))
      colnames(table1) <- c("Selected Files:")
      return(table1)
    } else if(!is.null(listObject$over.clustering)) {
      table1 <- as.matrix(names(listObject$over.clustering))
      colnames(table1) <- c("Selected Files:")
      return(table1)
    } else {
      return(NULL)
    }
    },colnames = T, width = "100%"
  )

  output$boxDelBut <- renderUI({
    if(is.null(listObject$flow.frames) && is.null(listObject$over.clustering)) return(NULL)
    box(title="Unselect Files", collapsible= TRUE, width = "100%",id="inputboxGroups", status="success",
      fluidRow(
        column(11, id = "align_button",tableOutput("flow_frames")),
        column(1,uiOutput("buttonDel"))
      )
    )
  })

  #Function creating ui buttons to remove files from the list of FCS files.
  observe({
    if(is.null(listObject$flow.frames) && is.null(listObject$over.clustering)){return(NULL)}
    if(is.null(listObject$flow.frames)){list <- listObject$over.clustering}
    if(is.null(listObject$over.clustering)){list <- listObject$flow.frames}
    listObject$files.id <- unlist(lapply(c(1:length(list)),function(x) {
        return(paste0(sample(letters, x + 1, replace = TRUE), collapse = ""))
    }))
    del_button_output <- lapply(c(1:length(listObject$files.id)), function(x) {
      del_button_name <- paste0("delButton_", listObject$files.id[x])
      del_button_object <-actionButton(del_button_name,"" ,icon = icon(name = "trash", lib = "glyphicon"))
      return(del_button_object)
    })
    do.call(tagList, del_button_output)
    output$buttonDel <- renderUI({
      del_button_output
    })
  })

  #Function triggering when a delete button is used, deleting the corresponding file from the list.
  observe({
    if(is.null(listObject$flow.frames) && is.null(listObject$over.clustering)){return(NULL)}
    lapply(listObject$files.id, function(i){
      observeEvent(input[[paste0("delButton_", i)]], {
        if (length(listObject$files.id) <= 1) {
          return(NULL)
        } 
        index <- which(i == listObject$files.id)
        if(is.null(listObject$flow.frames)){
          listObject$over.clustering <- listObject$over.clustering[-index]
        } else if(is.null(listObject$over.clustering)){
          listObject$flow.frames <- listObject$flow.frames[-index]
        }  
      })
    })
  })

  ################################################################################
  ## GROUP CLUSTERING

  output$clusteringui3 <- renderUI({
    if(is.null(listObject$flow.frames)) return(NULL)
    dd <- listObject$clustering.groups
    if(length(listObject$clustering.groups)==0){listObject$clustering.groups <- NULL}
    return(tagList(mapply(
      get_cluster_groups_table, dd, names(dd), SIMPLIFY = F
    )))
  })

  output$clusteringui2 <- renderUI({
    if(is.null(listObject$flow.frames)) return(NULL)
    selectInput("clusteringui_files_list",label = "File list",
      choices = names(listObject$flow.frames),selectize = F,
      multiple = T,width = "100%"
    )
  })

  observeEvent(input$clusteringui_add_all_groups, {
    for(name in names(listObject$flow.frames)){
      listObject$clustering.groups <- c(listObject$clustering.groups, setNames(list(name), name[1]))
    }
  })

  observe({
    key <- input$clusteringui_remove_clustering_group$key
    if(!is.null(key) && key != ""){
      isolate({
        listObject$clustering.groups[key] <- NULL
      })
    }
  })

  observeEvent(input$clusteringui_add_clustering_group, {
    if(is.null(input$clusteringui_files_list)) {return(NULL)}
    files_list <- isolate({
      input$clusteringui_files_list
    })
    listObject$clustering.groups <-
      c(listObject$clustering.groups, setNames(list(files_list), files_list[1]))
  })

  output$boxMultiGroup <- renderUI({
    if(is.null(listObject$flow.frames)) return(NULL)
    box(title="Unselect Files", collapsible= TRUE, width = "100%",id="inputboxGroups",status="primary",
      fluidRow(
        column(12,uiOutput("clusteringui2"),
          column(6,actionButton("clusteringui_add_clustering_group", "Add clustering group")),
          column(6,actionButton("clusteringui_add_all_groups", "Add all files into groups")),
          tags$br(),tags$br(),tags$br(),
          uiOutput("clusteringui3")
        )
      )
    )
  })

########################################################################################
############################ PRE PROCESS BOX AND UI ####################################

  output$output_process <- renderText("No Process")

  observe({
    if(!is.null(input$preprocess) && input$preprocess>=1){return(NULL)}
    if(is.null(listObject$flow.frames) || is.null(listObject$clustering.groups)){
      output$boxPreProcess <- NULL
      return(NULL)
    }
    if(!is.null(listObject$over.clustering)){
      output$boxPreProcess <- NULL
      return(NULL)
    } 
    if(input$tab == "clustering"){
      output$boxPreProcess <- renderUI({
        box(title="Pre Processing",collapsible = TRUE,status="warning", width = "100%",id="inputboxPreprocess",
          fluidRow(
            column(5,
              selectInput("method_transform","Choose a transformation",selected = "None",
                choices = c("None", "Asinh", "Logicle"))
            ),
            column(5,
              numericInput("args_transform","Co-factor(Asinh) / Biexp (Logicle)",min = NA, max = NA, step = NA, value=NULL)
            ),
            column(2,
              checkboxInput("comp","Compensate", value=TRUE)
            )
          ),
          fluidRow(
            column(5,
              actionButton("preprocess","Run Pre Process")
            ),
            column(5,
              verbatimTextOutput("output_process")
            )
          ),
          fluidRow(
            column(5,
              selectInput("marker_untrans","Select Marker Untransformed",choice = c(NULL,listObject$params.ori), multiple=TRUE)
            )
          )
        )
      })
    }
  })

  observeEvent(input$method_transform,{
    if(input$method_transform=="Asinh"){
      updateNumericInput(session, "args_transform", min = NA, max = NA, step = 1, value=5)
    } else if(input$method_transform=="Logicle") {
      updateNumericInput(session, "args_transform", min = NA, max = NA, step = NA, value = NULL)
    } else {
      updateNumericInput(session, "args_transform", min = NA, max = NA, step = NA, value = NULL)
    }
  })

  observeEvent(input$preprocess,{

    progress <- Progress$new()
    shinyjs::disable("inputboxPreprocess")
    progress$set(message="Preprocessing...", value = 1)
    
    args_transform <- input$args_transform
    if(is.na(input$args_transform)){args_transform <- NULL}
    
    flow.frames <- listObject$flow.frames
    listObject$flow.frames <- pre_process_fcs(flow.frames, 
      arg = args_transform, 
      transformation = input$method_transform, 
      compens = input$comp,
      marker_untrans = input$marker_untrans)
    listObject$doneComp <- TRUE
    output$output_process <- renderText("Processing Done !!!")
    progress$close()
  })

########################################################################################
############################ CLUSTERING BOX AND UI #####################################

  output$output_clustering <- renderText("No Clustering")

  observe({
    if(!is.null(input$run_clustering) && input$run_clustering>=1){return(NULL)}
    if(is.null(listObject$flow.frames) || is.null(listObject$clustering.groups)){
      output$boxRunClustering <- NULL
      return(NULL)
    }
    if(!is.null(listObject$over.clustering)){
      output$boxRunClustering <- NULL
      return(NULL)
    }
    if(!is.null(input$preprocess) && input$preprocess>=1){return(NULL)}
    output$boxRunClustering <- renderUI({
      box(title="Run Clustering and select Over Cluster",collapsible = TRUE, status = "info", width = "100%",id="inputBoxRunClustering",
        fluidRow(
          column(5,
            selectInput("marker_e","Select Enrichment Marker (and don't run clustering)",
              choice = c("",listObject$params.ori), selected = 1, multiple=F
            )
          ),
          column(5,
            actionButton("action_e","Compute Enrichment")
          )
        ),tags$hr(style="border: 1.5px solid DarkSlateBlue ;"),
        fluidRow(
          column(3,
            numericInput("clustering_k","Number of cluster",min = NA, max = NA, step = NA, value=200)
          ),
          column(3,
            numericInput("clustering_args","Sampling",min = NA, max = NA, step = NA, value=50)
          ),
          column(3,
            numericInput("ncores", "Number of Cores to use", 
              min = 1, max = (detectCores()-1), step = 1, value = (detectCores()-1)
            )
          )
        ),
        fluidRow(
          column(6,
            selectInput("clusteringui_file_for_markers","Load marker names from file",
              choices = c(names(listObject$flow.frames)),width = "100%"
            )
          ),
          column(6,
            selectInput("clusteringui_markers","Choose clustering marker",
              choices = c("", listObject$params.ori), multiple = T, width = "100%"
            )
          )
        ),
        fluidRow(
          column(3,
            actionButton("run_clustering","Run Clustering")
          ),
          column(3,
            verbatimTextOutput("output_clustering")
          ),
          column(6,
            actionButton("add_all_markers_clustering", "Add All")
          )  
        )
      )
    })
  })

  observe({
    if(!is.null(input$clusteringui_file_for_markers)&&grepl("*.fcs$", input$clusteringui_file_for_markers)){
      updateSelectInput(session,"clusteringui_markers", selected = "", choices = c("", listObject$params.ori))
    }
  })
  
  observeEvent(input$add_all_markers_clustering, {
    if(!is.null(input$clusteringui_file_for_markers)&&grepl("*.fcs$", input$clusteringui_file_for_markers)){
      updateSelectInput(session,"clusteringui_markers", selected = listObject$params.ori, choices = c("", listObject$params.ori))
    }
  })

  observeEvent(input$action_e, {
    if(length(listObject$flow.frames)!=length(listObject$clustering.groups)){
      showModal(modalDialog(title = "Important message", "Use this argument juste with same number of file and group",easyClose=TRUE))
      return(NULL)
    }
    if(input$marker_e == "" || is.null(input$marker_e)){
      return(NULL)
    } else if (!is.null(listObject$flow.frames)){
      shinyjs::disable("inputboxPreprocess")
      progress <- Progress$new()
      progress$set(message="Compute with Enrichment...", value=1)
      marker_e <- input$marker_e
      listObject$over.clustering <- lapply(c(1:length(listObject$flow.frames)), get_matrix_from_fcs, flow.frames = listObject$flow.frames, method = "mean", marker = marker_e, groups.clustering = NULL) #groups.clustering = listObject$clustering.groups)
      names(listObject$over.clustering) <- gsub(".fcs$",".txt",names(listObject$flow.frames))
      listObject$flow.frames <- NULL
      progress$close()
    }
  })
 
  observeEvent(input$run_clustering,{
    if(is.null(input$clusteringui_markers)) return(NULL)
    shinyjs::disable("inputboxPreprocess")
    progress <- Progress$new()
    progress$set(message="Clustering in progress...", value=1)

    args_transform <- input$args_transform
    if(is.na(input$args_transform)){args_transform <- 1}

    # print(lisgroups.clustering)

    result <- run_clustering(flow.frames = listObject$flow.frames, 
      methods = "CLARA",
      args = input$clustering_args,
      nb.cluster = input$clustering_k,
      params = input$clusteringui_markers,
      outputDir = listObject$outputDirectory,
      index = NULL,
      groups = listObject$clustering.groups,
      ncores = input$ncores,
      transComp = c(input$comp, input$method_transform, args_transform, listObject$doneComp),
      marker_untrans = input$marker_untrans
    )
    progress$close()
   
    listObject$origin.parameters <- colnames(listObject$flow.frames[[1]]) #Heatmap Idea
    listObject$flow.frames.enrich  <- result[[1]]
    listObject$flow.frames <- NULL
    listObject$over.clustering <- result[[2]]
    names(listObject$over.clustering) <- gsub(".fcs$",".txt", names(listObject$over.clustering))
    output$output_clustering <- renderText("Clustering Done !!")
  })
  
  output$ddlTabl <- downloadHandler(
    filename = function() {
      paste("output", "zip", sep=".")
    },
    content = function(fname){
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      for(i in c(1:length(listObject$over.clustering))) {
        name <- names(listObject$over.clustering)[i]
        path <- paste0(name,".csv")
        fs <- c(fs, path)
        write.csv(listObject$over.clustering[[i]], path,sep=",",row.names = FALSE)
      }
      zip(zipfile=fname, files=fs)
      setwd(listObject$inputDirectory)
    },
    contentType = "application/zip"
  )

  observeEvent(input$computeHeatmap, {
    if(is.null(listObject$over.clustering) || length(input$select.params.heatmap)<2){return(NULL)}
    
    output$heatmap <- renderPlotly({
      heatmaply(as.matrix(listObject$over.clustering[[input$select.over.clustering]][,input$select.params.heatmap]),
      colors = plasma(256, alpha = 1, begin = 0, end = 1, direction = 1),
      scale = "none", Colv = "row", dendrogram = "row",
      height=1200)
    })
  })
 
########################################################################################
################################ ANALYSE BOX AND UI INPUT ##############################

  output$boxGated <- renderUI({
    if(input$tab == "analysis"){
      box(title="Choose Gated FCS",collapsible = TRUE,status = "danger", solidHeader = TRUE, width = 250,
        fileInput("gatedInput", "",
                multiple = TRUE,
                accept = c(".fcs")
        ),
        actionButton("refreshGated","Refresh Gated")
      )
    }
  })

  observeEvent(input$refreshGated,{
    listObject$gated.flow.frames <- NULL
  })

  output$boxGatedTable <- renderUI({
    if(is.null(listObject$gated.flow.frames)) return(NULL)
    if(input$tab != "mapping"){
        files <- names(listObject$gated.flow.frames)
        box(title="List Gated Files",collapsible = TRUE, status = "danger", solidHeader = TRUE, width = 250,
        renderTable({
          files <- as.matrix(files)
          colnames(files) <- "Names"
          return(files)})
      )
    }
  })

  observeEvent(input$gatedInput,{
    progress <- Progress$new()
    progress$set(message = "Read Upload file", value = 1)
    filesname <- c(names(listObject$gated.flow.frames),as.vector(input$gatedInput$name))
    listObject$gated.datapaths <- input$gatedInput$datapath
    new.flow.frames <- lapply(as.vector(listObject$gated.datapaths), function(x) {
      return(read.FCS(x, emptyValue = FALSE))
    })
    listObject$gated.flow.frames  <- c(listObject$gated.flow.frames, new.flow.frames)
    names(listObject$gated.flow.frames) <- filesname
    progress$close()
  })

  observe({
    if(length(listObject$gated.flow.frames)>1){
      output$boxMapMarkers <- renderUI({
        box(width = "100%", solidHeader = TRUE, title = "Align markers", collapsible = TRUE,
            column(6,uiOutput("analysisui1"),
                   uiOutput("analysisuipanel1")
            ),
            column(6,uiOutput("analysisui2"),
                   uiOutput("analysisuipanel2")
            ),
            uiOutput("match_markers")
        )
      })
      output$analysisui1 <- renderUI({
        labels <- pData(parameters(listObject$gated.flow.frames[[1]]))[,2]
        params.names <- pData(parameters(listObject$gated.flow.frames[[1]]))[,1]
        labels[which(is.na(labels))] <- colnames(listObject$gated.flow.frames[[1]])[c(which(is.na(labels)))]
        names(params.names) <- labels
        selectInput("ref_marker","Landmark Markers", choices = params.names,multiple=TRUE)
      })
      output$analysisuipanel1 <- renderUI({
        wellPanel(returnOrder("mappingui_ref_markers_list", c("")), style = "background-color: #dd4b39")
      })
      output$analysisui2 <- renderUI({
        selectInput("ori_marker","Analysis Markers", choices=colnames(listObject$over.clustering[[1]]),multiple=TRUE)
      })
      output$analysisuipanel2 <- renderUI({
        wellPanel(returnOrder("mappingui_ori_markers_list", c("")), style = "background-color: #3c8dbc")
      })
      output$match_markers <- renderUI({
        actionButton("match_markers_button", "Match Landmark markers with Analysis markers", width = "100%")
      })
    }
  })
  
  observeEvent(input$match_markers_button, {
    labels <- pData(parameters(listObject$gated.flow.frames[[1]]))[,2]
    params.names <- pData(parameters(listObject$gated.flow.frames[[1]]))[,1]
    labels[which(is.na(labels))] <- colnames(listObject$gated.flow.frames[[1]])[c(which(is.na(labels)))]
    names(params.names) <- labels
    
    tab <- names(params.names[match(input$ref_marker, params.names)])
    tab2 <- colnames(listObject$over.clustering[[1]])
    add <- list()
    for (i in tab) {
      if (i %in% tab2) {
        add <- c(add, i)
      }
    }
    updateSelectInput(
      session,
      "ori_marker",
      selected = add,
      choices = c("", colnames(listObject$over.clustering[[1]]))
    )
  })
  
  observe({
    if(!is.null(input$preprocessgated) && input$preprocessgated>=1){return(NULL)}
    if(is.null(listObject$gated.flow.frames)){
      output$boxPreProcess <- NULL
      return(NULL)
    }
    if(input$tab == "analysis"){
      output$boxPreProcessGated <- renderUI({
        box(title="Pre Processing Gated Files",collapsible = TRUE, solidHeader = TRUE, width = "100%",id="inputboxPreprocessGated",
            fluidRow(
              column(5,
                     selectInput("method_transform_gated","Choose a transformation",selected = "None",
                                 choices = c("None", "Asinh", "Logicle"))
              ),
              column(5,
                     numericInput("args_transform_gated","Co-factor(Asinh) / Biexp (Logicle)",min = NA, max = NA, step = NA, value=NULL)
              ),
              column(2,
                     checkboxInput("comp_gated","Compensate", value=TRUE)
              )
            ),
            fluidRow(
              column(5,
                     actionButton("preprocessgated","Run Pre Process")
              ),
              column(5,
                     verbatimTextOutput("output_process_gated")
              )
            )
        )
      })
    }
  })
  
  observeEvent(input$preprocessgated,{
    
    progress <- Progress$new()
    progress$set(message="Preprocessing...", value = 1)
    shinyjs::disable("inputboxPreprocessGated")
    
    args_transform <- input$args_transform_gated
    if(is.na(input$args_transform_gated)){args_transform <- NULL}
    
    flow.frames <- listObject$gated.flow.frames
    listObject$gated.flow.frames <- pre_process_fcs(flow.frames, 
                                              arg = args_transform, 
                                              transformation = input$method_transform_gated, 
                                              compens = input$comp_gated, marker_untrans = NULL)
    output$output_process <- renderText("Processing Done !!!")
    progress$close()
  })

  observeEvent(input$ref_marker,{
    if(!is.null(input$ref_marker) && length(input$ref_marker > 0)) {
      
      labels <- pData(parameters(listObject$gated.flow.frames[[1]]))[,2]
      params.names <- pData(parameters(listObject$gated.flow.frames[[1]]))[,1]
      labels[which(is.na(labels))] <- colnames(listObject$gated.flow.frames[[1]])[c(which(is.na(labels)))]
      names(params.names) <- labels
      
      tab <- names(params.names[match(input$ref_marker, params.names)])
      updateReturnOrder(session,
        "mappingui_ref_markers_list",
        tab
      )
    }
  })

  observeEvent(input$ori_marker,{
    if(!is.null(input$ori_marker) && length(input$ori_marker > 0)) {
      updateReturnOrder(session,
        "mappingui_ori_markers_list",
        input$ori_marker
      )
    }
  })
  
  observe({
    if (!is.null(input$ref_marker) && !is.null(input$ori_marker) && (length(input$ref_marker) == length(input$ori_marker)) && length(input$ref_marker)>1) {
      output$boxParameters <- renderUI(box(solidHeader = TRUE, title = "Run Analysis", collapsible = TRUE, width = "100%",
        selectInput("map.clustedFiles.names","Select Reference Map",choices = names(listObject$over.clustering), selected=1, multiple=FALSE),
        checkboxInput(inputId = "inter_cluster", label = "Add inter-cluster connections", value = TRUE),
        actionButton("run_analysis", label = "Run Analysis")
      ))
    } else {
      output$boxParameters <- NULL
    }
  })
  
  observeEvent(input$run_analysis, {

    
    progress <- Progress$new()
    progress$set(message="Analysis...", value=1)
     
    result <- run_analysis_gated(
          listObject$gated.flow.frames,
          listObject$over.clustering,
          map.clustedFiles.names = input$map.clustedFiles.names,
          listObject$outputDirectory,
          listObject$clustering.groups,
          input$mappingui_ref_markers_list,
          input$mappingui_ori_markers_list,
          col.names.inter_cluster = NULL,
          ew_influence = NULL,
          inter_cluster.weight_factor = 0.7,
          inter.cluster.connections = input$inter_cluster,
          overlap_method = "repel"
    )
    listObject$scaffold.data <- result
    
    progress$close()
  })

########################################################################################
################################ MAPING EXPLORATIONS ###################################

  output$boxScaffoldMap <- renderUI({
    if(input$tab == "mapping"){
      box(title="Load Scaffold",collapsible = TRUE,status = "success", solidHeader = TRUE, width = 250,
        fileInput("mapInput","",multiple = FALSE,accept = c(".scaffold"))
      )
    } else {
      return(NULL)
    }
  })

  output$boxScaffoldMapUI <- renderUI({
    if(is.null(listObject$scaffold.data)) {return(NULL)}
    if(input$tab != "mapping") {
      return(NULL)
    } else {
      box(" ", collapsible=TRUE, width=12,
      tabBox(title = " ",id = "tabset1", width=12, 
        tabPanel(title = "Maps",
          selectizeInput("graphui_selected_graph","Choose a graph:", choices = c("")),
          selectizeInput("graphui_active_sample","Active sample",choices = c("All")),
          selectInput("graphui_marker","Nodes color:", choices = c("Default")),
          fluidRow(
            column(6,
              selectInput("graphui_stats_type","Stats type",choices = c("Ratio", "Difference"))
            ),
            column(6,
              selectInput("graphui_stats_relative_to","Stats relative to:",choices = c("Absolute"))
            )
          ),
          selectInput("graphui_color_scaling","Color scaling:",choices = c("global", "local")),
          fluidRow(
            column(6,
              selectInput("graphui_node_size","Nodes size:",choices = c("Proportional", "Default"))
            ),
            column(6,
              numericInput("graphui_min_node_size","Minimum node size",2,min = 0,max = 1000)
            )
          ),
          fluidRow(
            column(6,
              numericInput("graphui_max_node_size","Maximum node size",40,
                min = 0,max = 1000
              )
            ),
            column(6,
              numericInput("graphui_landmark_node_size", "Landmark node size",8,min = 0,max = 1000)
            )
          ),
          selectInput("graphui_display_edges","Display edges:",
              choices = c("All", "Highest scoring", "Inter cluster", "To landmark"))
        ),
        tabPanel(title="Colors",
          selectInput("graphui_color_number", "Number of colors", choices = c(2, 3)),
          fluidRow(
            column(6,colourpicker::colourInput("graphui_color_under", "Under:", value = "#FFFF00")
            ),
            column(6,colourpicker::colourInput("graphui_color_over", "Over:", value = "#0000FF")
            )
          ),
          fluidRow(
            column(4,colourpicker::colourInput("graphui_color_min", "Min:", value = "#E7E7E7")),
            column(4,
              conditionalPanel(
                condition = "input.graphui_color_number == 3",
                colourpicker::colourInput("graphui_color_mid", "Mid:", value = "#E7E7E7")
              )
            ),
            column(4,colourpicker::colourInput("graphui_color_max", "Max:", value = "#E71601"))
          ),
          conditionalPanel(
            condition = "input.graphui_color_number == 3",
            sliderInput("graphui_color_scale_mid","Color scale midpoint",
              min = 0.0,max = 5.0,value = 2.5,round = -2,step = 0.1,sep = ""
            )
          ),
          sliderInput( "graphui_color_scale_lim","Color scale limits",
                       min = 0.0,max = 5.0,value = c(0.0, 5.0),round = -2,step = 0.1,sep = ""
          ),
          fluidRow(
            column(6, numericInput("graphui_color_scale_min", "Color scale min:", 0)),
            column(6,numericInput("graphui_color_scale_max", "Color scale max:", 5))
          )
        ))
      )}
  })
  
  output$boxScaffoldMapDownload <- renderUI({
    if(is.null(listObject$scaffold.data)) {return(NULL)}
    if(input$tab == "mapping"){
      box(title = "Download events and freqs table", collapsible=TRUE, width=12,
        downloadButton(outputId = "CellsPerLandmarkDownload", label = "Download", icon("download")))
    } else {
      return(NULL)
    }
  })

  output$CellsPerLandmarkDownload <- downloadHandler(
    filename = function() {
      return(paste("outputTable", "csv", sep="."))
    },
    content = function(f){
      table <- get_cells_per_landmark_all_files(listObject$scaffold.data)
      write.csv(table, f, sep=",")
    },
    contentType = "application/csv"
  )
  
  observeEvent(input$mapInput,{
    progress <- Progress$new()
    progress$set(message = "Loading maps", value = 0.33)
    listObject$scaffold.data <- my_load(input$mapInput$datapath)
    progress$set(message = "Loading maps", value = 0.66)
    updateSelectInput(session,"graphui_selected_graph", choices = c("", names(listObject$scaffold.data$graphs)))
    progress$close()
  })

  observe({
    if(input$tab == "mapping"){
      updateSelectInput(session,"graphui_selected_graph", choices = c("", names(listObject$scaffold.data$graphs)))
    }
  })

  output$boxMap <- renderUI({
    if(is.null(listObject$scaffold.data)) return(NULL)
    box(title="SCAFFOLD Map",collapsible = TRUE,status="warning", width = "100%",id="",
      reactiveNetwork(outputId = "graphui_mainnet")
    )
  })
   
  output$boxTable <- renderUI({
    if(is.null(listObject$scaffold.data)) return(NULL)
    box(title = "Reference Table", collapsible = TRUE, collapsed = TRUE, status = "warning", width = "100%", id = "boxTable",
        dataTableOutput("graphui_table"))
  })
  
  output$graphui_table <- renderDataTable({
    if (!is.null(sc.data) &&
        !is.null(input$graphui_selected_graph) &&
        input$graphui_selected_graph != "")
    {
      if (is.null(input$graphui_selected_nodes) ||
          length(input$graphui_selected_nodes) == 0)
      {
        get_number_of_cells_per_landmark(listObject$scaffold.data,
                                         input$graphui_selected_graph)
      }
      else
      {
        get_summary_table(
          listObject$scaffold.data,
          input$graphui_selected_graph,
          input$graphui_selected_nodes
        )
      }
    }
  }, options = list(
    scrollX = "1000px",
    searching = FALSE,
    scrollY = "800px",
    paging = FALSE,
    info = FALSE,
    processing = FALSE
  ))

  output$graphui_mainnet <- reactive({
    ret <- get_main_graph()
    if (!is.null(ret))
    {
      ret$color <- get_color()
      ret$trans_to_apply <- isolate({
        input$graphui_cur_transform
      })
    }
    return(ret)
  })

  get_main_graph <- reactive({
    sc.data <- listObject$scaffold.data
    if (!is.null(sc.data) &&
        !is.null(input$graphui_selected_graph) &&
        input$graphui_selected_graph != "")
    {
      attrs <-
        get_numeric_vertex_attributes(sc.data, input$graphui_selected_graph)
      node.size.attr <-
        combine_marker_sample_name("popsize", input$graphui_active_sample)

      isolate({
        sel.marker <- NULL
        if (input$graphui_marker %in% attrs)
          sel.marker <- input$graphui_marker
        else
          sel.marker <- "Default"
        updateSelectInput(
          session,
          "graphui_marker",
          choices = c("Default", attrs),
          selected = sel.marker
        )
        updateSelectInput(
          session,
          "graphui_markers_to_plot",
          choices = attrs,
          selected = attrs
        )
        sample.names <-
          get_sample_names(sc.data, input$graphui_selected_graph)
        updateSelectInput(
          session,
          "graphui_active_sample",
          choices = c("All", sample.names),
          selected = input$graphui_active_sample
        )
        updateSelectInput(
          session,
          "graphui_stats_relative_to",
          choices = c("Absolute", sample.names),
          selected = input$graphui_stats_relative_to
        )
      })
      return(
        get_graph(
          sc.data,
          input$graphui_selected_graph,
          node.size.attr,
          input$graphui_min_node_size,
          input$graphui_max_node_size,
          input$graphui_landmark_node_size
        )
      )
    }
    else
      return(NULL)
  })

  read_color_scale_info <- reactive({
    return(
      list(
        sel.marker = input$graphui_marker,
        color.scale.lim = input$graphui_color_scale_lim,
        color.scale.mid = input$graphui_color_scale_mid
      )
    )
  })

  get_color_scale <- reactive({
    #This code only updates the color scales
    sc.data <- listObject$scaffold.data
    if (is.null(sc.data) || is.null(get_main_graph()))
      return(NULL)
    sel.marker <- input$graphui_marker
    rel.to <- input$graphui_stats_relative_to
    color.scaling <- input$graphui_color_scaling
    stats.type <- input$graphui_stats_type
    isolate({
      color <- NULL
      if (sel.marker != "")
      {
        #Colors are not really important here, only included because they need to be passed to the function
        min.color <- input$graphui_color_min
        mid.color <- input$graphui_color_mid
        max.color <- input$graphui_color_max
        under.color <- input$graphui_color_under
        over.color <- input$graphui_color_over
        color <-
          get_color_for_marker(
            sc.data,
            sel.marker,
            rel.to,
            input$graphui_selected_graph,
            input$graphui_active_sample,
            color.scaling,
            stats.type,
            colors.to.interpolate = c(min.color, mid.color, max.color),
            under.color,
            over.color
          )
        if (!is.null(color$color.scale.lim)
            && !(is.null(color.scaling)) && color.scaling == "local")
        {
          updateSliderInput(
            session,
            "graphui_color_scale_lim",
            min = color$color.scale.lim$min,
            max = color$color.scale.lim$max,
            step = 0.1,
            value = c(
              color$color.scale.lim$min,
              color$color.scale.lim$max
            )
          )
          updateSliderInput(
            session,
            "graphui_color_scale_mid",
            min = color$color.scale.lim$min,
            max = color$color.scale.lim$max,
            step = 0.1,
            value = mean(
              c(
                color$color.scale.lim$min,
                color$color.scale.lim$max
              )
            )
          )
        }
      }
    })
  })

  get_color <- reactive({
    #This code does the actual coloring
    get_color_scale()
    color.scale.info <- read_color_scale_info()
    min.color <- input$graphui_color_min
    mid.color <- input$graphui_color_mid
    max.color <- input$graphui_color_max
    under.color <- input$graphui_color_under
    over.color <- input$graphui_color_over
    color.scale.lim <- color.scale.info$color.scale.lim
    colors.to.interpolate <- NULL
    color.scale.mid <- NULL
    if (input$graphui_color_number == 3)
    {
      colors.to.interpolate <- c(min.color, mid.color, max.color)
      color.scale.mid <- color.scale.info$color.scale.mid
    }
    else
      colors.to.interpolate <- c(min.color, max.color)
    return(isolate({
      sel.marker <- color.scale.info$sel.marker

      color.vector <- NULL
      active.sample <- input$graphui_active_sample
      rel.to <- input$graphui_stats_relative_to
      color.scaling <- input$graphui_color_scaling
      stats.type <- input$graphui_stats_type

      if (sel.marker != "")
      {
        sc.data <- listObject$scaffold.data
        if (!is.null(sc.data))
        {
          color <-
            get_color_for_marker(
              sc.data,
              sel.marker,
              rel.to,
              input$graphui_selected_graph,
              active.sample,
              color.scaling,
              stats.type,
              colors.to.interpolate = colors.to.interpolate,
              under.color,
              over.color,
              color.scale.limits = color.scale.lim,
              color.scale.mid = color.scale.mid
            )
          color.vector <- color$color.vector
        }
      }
      return(color.vector)
    }))
  })

  output$graphui_table <- renderDataTable({
    if (!is.null(listObject$scaffold.data) &&
        !is.null(input$graphui_selected_graph) &&
        input$graphui_selected_graph != "")
    {
      if (is.null(input$graphui_selected_nodes) ||
          length(input$graphui_selected_nodes) == 0)
      {
        get_number_of_cells_per_landmark(listObject$scaffold.data,
          input$graphui_selected_graph)
      }
      else
      {
        get_summary_table(
          listObject$scaffold.data,
          input$graphui_selected_graph,
          input$graphui_selected_nodes
        )
      }
    }
    },options = list(scrollX = TRUE,searching = FALSE,scrollY = "800px",
      paging = FALSE,
      info = FALSE,
      processing = FALSE
  ))

  output$graphui_dialog1 <- reactive({
    sc.data <- listObject$scaffold.data
    ret <- ""
    if (!is.null(sc.data))
      ret <-
      sprintf("Markers used for SCAFFoLD: %s",
              paste(sc.data$scaffold.col.names, collapse = ", "))
    return(ret)
  })

  output$graphui_plot <- renderPlot({
    p <- NULL
    if (!is.null(input$graphui_plot_clusters) &&
        input$graphui_plot_clusters != 0)
    {
      isolate({
        col.names <- input$graphui_markers_to_plot
        if ((length(col.names) >= 1) &&
            (length(input$graphui_selected_nodes) >= 1))
          p <-
            plot_cluster(
              listObject$scaffold.data,
              input$graphui_selected_nodes,
              input$graphui_selected_graph,
              input$graphui_markers_to_plot,
              input$graphui_pool_cluster_data,
              input$graphui_plot_type
            )
      })
    }
  })

  observe({
    if (!is.null(input$graphui_reset_colors) &&
        input$graphui_reset_colors != 0)
    {
      session$sendCustomMessage(type = "reset_colors", "none")
    }
  })

  observe({
    if (!is.null(input$graphui_reset_graph_position) &&
        input$graphui_reset_graph_position != 0)
    {
      session$sendCustomMessage(type = "reset_graph_position", "none")
    }
  })

  observe({
    if (!is.null(input$graphui_toggle_landmark_labels) &&
        input$graphui_toggle_landmark_labels != 0)
    {
      display <-
        ifelse(input$graphui_toggle_landmark_labels %% 2 == 0,
               "",
               "none")
      session$sendCustomMessage(type = "toggle_label", list(target = "landmark", display = display))
    }
  })

  observe({
    display_edges <- input$graphui_display_edges
    session$sendCustomMessage(type = "toggle_display_edges", display_edges)
  })

  observe({
    if (!is.null(input$graphui_toggle_cluster_labels) &&
        input$graphui_toggle_cluster_labels != 0)
    {
      display <-
        ifelse(input$graphui_toggle_cluster_labels %% 2 == 0,
               "none",
               "")
      session$sendCustomMessage(type = "toggle_label", list(target = "cluster", display = display))
    }
  })

  observe({
    display <- tolower(input$graphui_node_size)
    session$sendCustomMessage(type = "toggle_node_size", list(display = display))
  })

  observe({
    if (!is.null(input$graphui_toggle_node_size) &&
        input$graphui_toggle_node_size != 0)
    {
      display <-
        ifelse(input$graphui_toggle_node_size %% 2 == 0,
               "proportional",
               "default")
      session$sendCustomMessage(type = "toggle_node_size", list(display = display))
    }
  })

########################################################################################
################################### MAPING DATASETS ####################################

  output$boxScaffoldMapData <- renderUI({
    if(input$tab == "mapdata" || input$tab == "export"){
      box(title="Load Scaffold Map",collapsible = TRUE,status = "success", solidHeader = TRUE, width = 250,
          fileInput("mapDataInput","",multiple = FALSE,accept = c(".scaffold"))
      )
    }
  })

  output$boxScaffoldMapView <- renderUI({
    if(!is.null(listObject$scaffold.data)){
      tmp <- names(listObject$scaffold.data$graphs)
      tmp <- as.matrix(tmp)
      colnames(tmp) <- "graphs"
      box(title="Map in your scaffold file",collapsible = TRUE,status = "success", solidHeader = TRUE, width = 250,
        renderTable({tmp})
      )
    }
  })

  observeEvent(input$mapDataInput, {
    progress <- Progress$new()
    progress$set(message = "Loading maps", value = 0.33)
    listObject$scaffold.data <- my_load(input$mapDataInput$datapath)
    progress$set(message = "Loading maps", value = 0.66)
    updateSelectInput(session,"graphui_selected_graph", choices = c("", names(listObject$scaffold.data$graphs)))
    progress$close()
  })
  
  observe({
    if(!is.null(listObject$scaffold.data)){
      output$boxMapDataMarkers <- renderUI({
        box(width = "100%", solidHeader = TRUE, title = "Align markers", collapsible = TRUE,
            column(6,uiOutput("mapdataui1"),
                   uiOutput("mapdatauipanel1")
            ),
            column(6,uiOutput("mapdataui2"),
                   uiOutput("mapdatauipanel2")
            ),
            uiOutput("match_markers_data")
        )
      })
      output$mapdataui1 <- renderUI({
        selectInput("ref_marker_data","Map Markers", choices=listObject$scaffold.data$scaffold.col.names,multiple=TRUE)
      })
      output$mapdatauipanel1 <- renderUI({
        wellPanel(returnOrder("mapdataui_ref_markers_list", NULL), style = "background-color: #00a65a")
      })
      output$mapdataui2 <- renderUI({
        selectInput("ori_marker_data","Analysis Markers", choices=colnames(listObject$over.clustering[[1]]),multiple=TRUE)
      })
      output$mapdatauipanel2 <- renderUI({
        wellPanel(returnOrder("mapdataui_ori_markers_list", NULL), style = "background-color: #3c8dbc")
      })
      output$match_markers_data <- renderUI({
        actionButton("match_markers_data_button", "Add All Map markers and match with Clustered markers", width = "100%")
      })
    }
  })
  
  observeEvent(input$match_markers_data_button, {
    tab <- listObject$scaffold.data$scaffold.col.names
    tab2 <- colnames(listObject$over.clustering[[1]])
    add <- list()
    for (i in tab) {
      if (i %in% tab2) {
        add <- c(add, i)
      }
    }
    updateSelectInput(
      session,
      "ref_marker_data",
      selected = tab,
      choices = c("", tab)
    )
    updateSelectInput(
      session,
      "ori_marker_data",
      selected = add,
      choices = c("", colnames(listObject$over.clustering[[1]]))
    )
  })
  
  observeEvent(input$ref_marker_data,{
    if(!is.null(input$ref_marker_data) && length(input$ref_marker_data)>0) {
      updateReturnOrder(session,
                        "mapdataui_ref_markers_list",
                        input$ref_marker_data
      )
    }
  })
  
  observeEvent(input$ori_marker_data,{
    if(!is.null(input$ori_marker_data) && length(input$ori_marker_data)>0) {
      updateReturnOrder(session,
                        "mapdataui_ori_markers_list",
                        input$ori_marker_data
      )
    }
  })
  
  observe({
    if (!is.null(input$ref_marker_data) && !is.null(input$ori_marker_data) && (length(input$ref_marker_data) == length(input$ori_marker_data)) && length(input$ref_marker_data)>1) {
      output$boxMapDataParams <- renderUI(box(solidHeader = TRUE, title = "Run Mapping", collapsible = TRUE, width = "100%",
                                           selectInput(inputId = "mapping_method", label = "Mapping Method", selected = "Concatenation", choices = c("Concatenation", "From Scratch")),
                                           checkboxInput(inputId = "inter_cluster_data", label = "Add inter-cluster connections", value = TRUE),
                                           actionButton("run_mapping", label = "Run Mapping")
      ))
    } else {
      output$boxParameters <- NULL
    }
  })
  
  #Runs the analysis and processes the files.
  observeEvent(input$run_mapping, {
    progress <- Progress$new()
    progress$set(message="Mapping...", value=1)
    
    if (!is.null(listObject$flow.frames)) 
      flow.frames <- listObject$flow.frames
    else if (!is.null(listObject$flow.frames.e))
      flow.frames <- listObject$flow.frames.e
    else
      flow.frames <- NULL
    
    col.names <- input$mapdataui_ori_markers_list
    ref.col.names <- input$mapdataui_ref_markers_list
    names.map <- ref.col.names
    #Missing values (i.e. non-mapped markers) are filled with NA
    names(names.map) <- col.names
    ew_influence <- NULL
    
    result <- run_analysis_existing(
      listObject$scaffold.data,
      names(listObject$over.clustering)[[1]],
      listObject$over.clustering,
      listObject$outputDirectory,
      col.names.matrix = input$mapdataui_ori_markers_list,
      col.names.map = input$mapdataui_ref_markers_list,
      inter.cluster.connections = input$mappingui_inter_cluster_connections,
      mode = input$mapping_method,
      names.map = names.map,
      col.names.inter_cluster = NULL,
      inter_cluster.weight_factor = NULL,
      overlap_method = "repel",
      ew_influence = ew_influence
      
    )
    print("Mapping done.")
    listObject$scaffold.data <- result
    progress$close()
  })

########################################################################################
################################### EXPORT DATASETS ####################################
  
  output$boxScaffoldMapExport <- renderUI({
    if(input$tab == "export"){
      box(title="Load Scaffold Map",collapsible = TRUE,status = "success", solidHeader = TRUE, width = 250,
          fileInput("mapExpInput","",multiple = FALSE,accept = c(".scaffold"))
      )
    }
  })
  
  observeEvent(input$mapExpInput, {
    progress <- Progress$new()
    progress$set(message = "Loading maps", value = 0.33)
    listObject$scaffold.data <- my_load(input$mapExpInput$datapath)
    progress$set(message = "Loading maps", value = 0.66)
    updateSelectInput(session,"graphui_selected_graph", choices = c("", names(listObject$scaffold.data$graphs)))
    progress$close()
  })
  
  output$boxScaffoldMapPops <- renderUI({
    if(!is.null(listObject$scaffold.data))
      if(input$tab == "export"){
        table.node <- scaffold_node_export(listObject$scaffold.data)
        box(title="Population's List", collapsible=TRUE, status = "success",solidHeader = TRUE, width=250,
            div(style = 'overflow-x:scroll;',renderTable({
              return(table.node)
            })),
            tags$br(),
            downloadButton("scaffoldNode","Pop Scaffold")
        )
      }
  })
  
  output$scaffoldNode <- downloadHandler(
    filename = function(){
      return("scaffoldNode.csv")
    },
    content = function(filename){
      data <- scaffold_node_export(listObject$scaffold.data)
      write.csv(data, filename)
    }
  )  
  
  observe({
    if(is.null(listObject$over.clustering) && is.null(listObject$flow.frames.enrich) && is.null(listObject$flow.frames)) return(NULL)
    if(!is.null(listObject$scaffold.data)){
      if(is.null(listObject$over.clustering)){data <- listObject$flow.frames}
      if(is.null(listObject$flow.frames)){data <- listObject$over.clustering}
            
      output$boxExportNames <- renderUI({
        box(title="Select names Map and Files", collapsible = TRUE, solidHeader=FALSE, status = "success", width="100%",
            fluidRow(
              column(6,
                     selectInput("map_files","Select Map Files", choices = names(listObject$scaffold.data$graphs), multiple=TRUE),
                     wellPanel(returnOrder("mapdataui_map_list", NULL), style = "background-color: #00a65a")
              ),
              column(6,
                     selectInput("names_files","Select Files", choices = names(data),multiple=TRUE),
                     wellPanel(returnOrder("mapdataui_files_list", NULL), style = "background-color: #3c8dbc")
              )
            ),
            fluidRow(
              column(6,
                actionButton("add_all_maps","Add All Maps")
              ),
              column(6,
                actionButton("add_all_files", "Add All Files")
              )
            )
        )
      })
    }
  })
  
  observeEvent(input$add_all_maps,{
    updateSelectInput(session, "map_files", selected = names(listObject$scaffold.data$graphs))
  })

  observeEvent(input$add_all_files,{
    if(is.null(listObject$over.clustering)){data <- listObject$flow.frames}
    if(is.null(listObject$flow.frames)){data <- listObject$over.clustering}
    updateSelectInput(session, "names_files", selected = names(data))
  })

  observeEvent(input$map_files,{
    if(!is.null(input$map_files)) {
      updateReturnOrder(session,"mapdataui_map_list",input$map_files)
    }
  })
  
  observeEvent(input$names_files,{
    if(!is.null(input$names_files)) {
      updateReturnOrder(session,"mapdataui_files_list", input$names_files)
    }
  })
  
  output$boxSelectEnrichment <- renderUI({
    if(is.null(listObject$scaffold.data)) return(NULL)
    if(length(input$mapdataui_map_list) != length(input$mapdataui_files_list)) {return(NULL)}
    if(length(input$mapdataui_map_list)<1 || length(input$mapdataui_files_list)<1) {return(NULL)}
    if(is.null(listObject$over.clustering) && is.null(listObject$flow.frames.enrich) && is.null(listObject$flow.frames)) return(NULL)
    if(input$tab == "export"){
      if(is.null(listObject$over.clustering)){data <- listObject$flow.frames}
      if(is.null(listObject$flow.frames)){data <- listObject$over.clustering}
      box(title="Export Scaffold MAP", collapsible = TRUE, solidHeader=FALSE, status = "success", width="100%",
          fluidRow(
            column(8,
                   selectInput("clusterID","Select Params Annotation (use just with enrichment fcs)",choices=colnames(data[[1]]), multiple=FALSE) 
            ),
            column(4,downloadButton("ddlAnot", "DDL Anotation"))
          ),
          fluidRow(
            column(4,selectInput("methodsMFI","MFIs Method",choices=c("median","mean"))),
            column(4,selectInput("sourceMFI","MFIs Compute for",choices=c("cluster","populations"))),
            column(4,downloadButton("ddlMFI","Download MFI"))
          )
      )
    } else {
      return(NULL)
    }
  })
  
  output$ddlAnot <- downloadHandler(
    filename = function(){
      return("output.zip")
    },
    content = function(file){
      progress <- Progress$new()
      progress$set(message="Annotation progress",value=1)
      
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      
      if(is.null(listObject$flow.frames.enrich) && is.null(listObject$flow.frames)){
        tables <- scaffold_cluster_export(
          input$mapdataui_map_list,
          input$mapdataui_files_list,
          listObject$over.clustering, 
          listObject$scaffold.data
        )
        for(i in c(1:length(tables))){
          name <- names(listObject$over.clustering)[i]
          path <- gsub(".csv$","_annotCluster.csv",name)
          fs <- c(fs, path)
          write.csv(tables[[i]], path, sep=",", row.names = FALSE)
        }
        
      } else if(is.null(listObject$over.clustering)){
        flow.frames <- scaffold_events_export(
          input$mapdataui_map_list,
          input$mapdataui_files_list,
          listObject$flow.frames, 
          listObject$scaffold.data, 
          input$clusterID
        )
        for(i in c(1:length(flow.frames))){
          name <- names(listObject$flow.frames)[i]
          fcs <- flow.frames[[i]]
          path <- gsub(".fcs$","_cellType.fcs",name)
          fs <- c(fs, path)
          write.FCS(fcs, path, delimiter="#")
        }
      }
      
      zip(zipfile=file, files=fs)
      setwd(listObject$inputDirectory)
      progress$close()
    },
    contentType = "application/zip"
  )
  
  output$ddlMFI <- downloadHandler(
    filename=function(){
      return("output.zip")
    },
    content = function(file){
      progress <- Progress$new()
      progress$set(message="Compute MFI",value=1)
      
      fs <- c()
      tmpdir <- tempdir()
      setwd(tempdir())
      
      if(is.null(listObject$flow.frames.enrich) && is.null(listObject$flow.frames)){
        # tables <- scaffold_cluster_export(
        #   input$mapdataui_map_list,
        #   input$mapdataui_files_list,
        #   listObject$over.clustering, 
        #   listObject$scaffold.data
        # )
        # for(i in c(1:length(tables))){
        #   name <- names(listObject$over.clustering)[i]
        #   path <- gsub(".csv$","_annotCluster.csv",name)
        #   fs <- c(fs, path)
        #   write.csv(tables[[i]], path, sep=",", row.names = FALSE)
        # }
        
      } else if(is.null(listObject$over.clustering)){
        
        print(listObject$flow.frames)
        print(input$mapdataui_map_list)
        print(input$mapdataui_files_list)
        mfi.table <- scaffold_pop_mfi(
          input$mapdataui_map_list,
          input$mapdataui_files_list,
          listObject$flow.frames, 
          listObject$scaffold.data, 
          input$clusterID,
          input$methodsMFI
        )
        
        print(length(mfi.table))
        for(i in c(1:length(mfi.table))){
          name <- names(listObject$flow.frames)[i]
          mat <-  mfi.table[[i]]
          path <- gsub(".fcs$","_MFIpop.csv",name)
          fs <- c(fs, path)
          write.csv(mat, path, row.names = FALSE)
        }
      }
      
      zip(zipfile=file, files=fs)
      setwd(listObject$inputDirectory)
      progress$close()
    },
    contentType = "application/zip"
  )
  
  
  ############################## NEWS ##############################
  
  # observe({
  #   if(length(listObject$flow.frames) < 1){return(NULL)}
  #   listObject$files.id <- as.vector(unlist(lapply(c(1:length(listObject$flow.frames)), function(x) {
  #       return(paste0(sample(letters, x + 1, replace = TRUE), collapse = ""))
  #    })))
  
  #   del_button_output <- lapply(c(1:length(listObject$files.id)), function(x) {
  #   	del_button_name <- paste0("delButton_", listObject$files.id[x])
  #   	del_button_object <- actionButton(del_button_name,"",icon = icon(name = "trash", lib = "glyphicon"))
  #     return(del_button_object)
  #   })
  
  #   do.call(tagList, del_button_output)
  #   output$buttonDel <- renderUI({
  #     del_button_output
  #   })
  # })
  
  
  # shinyFileChoose(input, "fcs_file", roots=roots, filetypes=c("fcs"))

})
