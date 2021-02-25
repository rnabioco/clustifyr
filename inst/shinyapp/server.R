# Define server logic to read selected file ----
server <- function(input, output, session) {
  # parse url to direct tab
  observe({
    query <- parseQueryString(session$clientData$url_search)
  
    if (!is.null(query[['tab']])) {
      updateTabItems(session, "tabs", query[['tab']])
    }
  })
  
  # reactive file location to make interactivity easier
  rv <- reactiveValues()
  rv$matrixloc <- NULL
  rv$metaloc <- NULL
  rv$step <- 0
  rv$clustifym <- "clustifyr not yet run"
  rv$lastgeo <- "GSE151974"#GSE113049"
  rv$ref <- NULL
  rv$ref_visited <- 0
  rv$ref_link <- NULL
  rv$res_visited <- 0
  rv$obj <- NULL
  rv$links2 <- data.frame()

  # waiter checkpoints
  w1 <- Waiter$new(
    id = "contents1",
    html = tagList(
      spin_flower(),
      h4("Matrix loading..."),
      h4("")
    )
  )
  
  w2 <- Waiter$new(
    id = "contents2",
    html = tagList(
      spin_flower(),
      h4("Metadata loading..."),
      h4("")
    )
  )
  
  w3 <- Waiter$new(
    id = "reference",
    html = tagList(
      spin_flower(),
      h4("Reference building..."),
      h4("")
    )
  )
  
  w4 <- Waiter$new(
    id = "clustify",
    html = tagList(
      spin_flower(),
      h4("Clustifyr running..."),
      h4("")
    )
  )
  
  w5 <- Waiter$new(
    id = "hmap",
    html = tagList(
      spin_flower(),
      h4("Heatmap drawing..."),
      h4("")
    )
  )

  w6 <- Waiter$new(id = "modalgeo",
                   html = tagList(
                     spin_flower(),
                     h4("Info fetching..."),
                     h4("")
                   ))

  w7 <- Waiter$new(id = "modalfiles",
                   html = tagList(
                     spin_flower(),
                     h4("File previewing..."),
                     h4("")
                   ))

  w8 <- Waiter$new(
    id = "contents3",
    html = tagList(
      spin_flower(),
      h4("Reference loading..."),
      h4("")
    )
  )

  data1 <- reactive({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.

    if (!is.null(input$file1) | !is.null(rv$matrixloc)) {
      if (!is.null(input$file1)) {
        rv$matrixloc <- input$file1
      }
      file <- rv$matrixloc

      if (!is.null(file)) {
        w1$show()
        message(file)
      }
      
      fileTypeFile1 <- tools::file_ext(file$datapath)
      req(file)
      if (str_to_lower(fileTypeFile1) == "rds") {
        df1 <- readRDS(file$datapath) 
        if (any(class(df1) %in% c("SingleCellExperiment", "Seurat"))) {
          rv$obj <- df1
          df1 <- object_data(rv$obj, "data")
        }
      } else if (str_to_lower(fileTypeFile1) == "rdata" | str_to_lower(fileTypeFile1) == "rda") {
        df1 <- load_rdata(file$datapath)
        if (any(class(df1) %in% c("SingleCellExperiment", "Seurat"))) {
          rv$obj <- df1
          df1 <- object_data(df1, "data")
        }
      } else {
        df1 <- fread(file$datapath)
      }
    } else if (!is.null(rv$obj)) {
      df1 <- object_data(rv$obj, "data")
    } else {
      return(NULL)
    }
    
    if ((!is_local) & (object.size(df1) > 3e9)) {
      message("Potential memory issue due to file size")
      showModal(modalDialog(
        h2("File size over 3Gb, online version may be unstable"),
        easyClose = TRUE,
        fade = FALSE,
        footer = NULL
      ))
    }
    
    df1 <- df1 %>% as.data.frame()
    if (!has_rownames(df1) & length(unique(df1[, 1])) == nrow(df1)) {
      rownames(df1) <- df1[, 1]
      df1[, 1] <- NULL
    }
    
    w1$hide()
    df1
  })
  
  data2 <- reactive({
    if (!is.null(input$file2) | !is.null(rv$metaloc)) {
      if (!is.null(input$file2)) {
        rv$metaloc <- input$file2
      }
      file <- rv$metaloc
      
      if (!is.null(file)) {
        w2$show()
        message(file)
      }
      
      fileTypeFile2 <- tools::file_ext(file$datapath)
      req(file)
      if (str_to_lower(fileTypeFile2) == "rds") {
        df2 <- readRDS(file$datapath)
        if (any(class(df2) %in% c("SingleCellExperiment", "Seurat"))) {
          rv$obj <- df2
          df2 <- object_data(df2, "meta.data")
        }
      } else if (str_to_lower(fileTypeFile2) == "rdata" | str_to_lower(fileTypeFile2) == "rda") {
        df2 <- load_rdata(file$datapath)
        if (any(class(df2) %in% c("SingleCellExperiment", "Seurat"))) {
          rv$obj <- df2
          df2 <- object_data(rv$obj, "meta.data")
        }
      } else {
        df2 <- fread(file$datapath)
      }
    } else if (!is.null(rv$obj)) {
      df2 <- object_data(rv$obj, "meta.data")
    } else {
      return(NULL)
    }

    if ((!is_local) & (object.size(df2) > 3e9)) {
      message("Potential memory issue due to file size")
      showModal(modalDialog(
        h2("File size over 3Gb, online version may be unstable"),
        easyClose = TRUE,
        fade = FALSE,
        footer = NULL
      ))
    }
    
    df2 <- df2 %>% as.data.frame()
    if (!has_rownames(df2) & length(unique(df2[, 1])) == nrow(df2)) {
      rownames(df2) <- df2[, 1]
      df2[, 1] <- NULL
    }

    w2$hide()
    df2
  })
  
  data3a <- reactive({
    if (!is.null(input$file3)) {
      rv$ref <- input$file3
    }
    if (rv$ref == "built-in") {
      return(NULL)
    }
    file <- rv$ref

    if (!is.null(file)) {
      w8$show()
      message(file)
    }

    fileTypeFile3 <- tools::file_ext(file$datapath)
    req(file)

    if (str_to_lower(fileTypeFile3) == "rds") {
      df3 <- readRDS(file$datapath) %>% as.data.frame()
    } else if (str_to_lower(fileTypeFile3) == "rdata" | str_to_lower(fileTypeFile3) == "rda") {
      df3 <- load_rdata(file$datapath) %>% as.data.frame()
    } else {
      df3 <- fread(file$datapath) %>% # , header = input$header, sep = input$sepMat) %>%
        as.data.frame()
    }

    if ((!is_local) & (object.size(df3) > 3e9)) {
      message("Potential memory issue due to file size")
      showModal(modalDialog(
        h2("File size over 3Gb, online version may be unstable"),
        easyClose = TRUE,
        fade = FALSE,
        footer = NULL
      ))
    }
    
    if (!has_rownames(df3) & length(unique(df3[, 1])) == nrow(df3)) {
      rownames(df3) <- df3[, 1]
      df3[, 1] <- NULL
    }

    w8$hide()
    df3
  })

  output$contents1 <- DT::renderDataTable({
    if (is.null(rv$matrixloc) & is.null(rv$obj)) {
      return(df1 <- data.frame(`nodata` = rep("", 6)))
    } else {
      df1 <- data1()
    }

    # file 1
    if (input$dispMat == "head") {
      cols <- ncol(df1)
      df1 <- df1[, 1:min(cols, 5)]
      return(head(df1))
    }
    else {
      return(df1)
    }
  })
  
  output$contents2 <- DT::renderDataTable({
    if (is.null(rv$metaloc) & is.null(rv$obj)) {
      return(df2 <- data.frame(`nodata` = rep("", 6)))
    } else {
      df2 <- data2()
    }
    
    updateSelectInput(session, "metadataCellType",
                      choices = c("", colnames(df2)),
                      selected = ""
    )

    # file 2
    if (input$dispMeta == "head") {
      return(head(df2))
    }
    else {
      return(df2)
    }

  },
  callback = DT::JS(js),
  selection = list(target = 'column', mode = "single"))

  output$colclicked <- renderUI({
    if (is.null(input[["column_clicked"]])) {
      "please select cluster column in drop-down menu, or click in the table"
    } else {
      input$metadataCellType
    }
  })

  observeEvent(input[["column_clicked"]], {
    updateSelectInput(session, "metadataCellType", 
                      selected = input[["column_clicked"]]                   
    )
  })

  output$ref_summary <- renderUI({
    HTML(paste0("cell types: ", ncol(data3()),
                "<br>",
                "genes: ", nrow(data3())))
  })

  data3b <- reactive({
    w8$show()
    rv$ref <- "built-in"
    ref <- refs[[ref_dict[input$dataHubReference]]]
    rv$ref_link <- refs_meta[ref_dict[input$dataHubReference], ] %>% pull(sourceurl)
    w8$show()

    ref
  })

  data3 <- reactive({
    b <- data3b()
    a <- data3a()
    if (is.null(rv$ref)) {
      return(data.frame(`nodata` = rep("", 6)))
    } else if (rv$ref == "built-in") {
      df3 <- b
    } else {
      df3 <- a
    }
    df3
  })

  output$contents3 <- DT::renderDataTable({
    df3 <- data3()

    # file 3
    if (input$dispMat == "head") {
      return(head(df3))
    }
    else {
      return(df3)
    }
  })
  
  observeEvent(input$matrixPopup, {
    showModal(modalDialog(
      tags$caption("Matrix table"),
      DT::renderDataTable({
        matrixRender <- head(data1())
        DT::datatable(matrixRender)
      }),
      easyClose = TRUE
    ))
  })
  
  observeEvent(input$metadataPopup, {
    showModal(modalDialog(
      tags$caption("Metadata table"),
      DT::renderDataTable({
        matrixRender <- head(data2())
        DT::datatable(matrixRender)
      }),
      easyClose = TRUE
    ))
  })

  data_avg <- reactive({
    if (input$metadataCellType == "") {
      return(NULL)
    }
    w3$show()
    reference_matrix <- average_clusters(mat = data1(), metadata = data2()[[input$metadataCellType]], if_log = FALSE)
    w3$hide()
    reference_matrix
  })
  
  dataClustify <- reactive({
    if (input$metadataCellType == "") {
      return(NULL)
    }
    w4$show()
    benchmarkRef <- data3()
    
    if (!is.null(rv$obj)) {
      message("Single cell object detected")
      matrixSeuratObject <- rv$obj
      if (any(class(matrixSeuratObject) == "SingleCellExperiment")) {
        matrixSeuratObject <- as.Seurat(matrixSeuratObject)
      }
    } else {
      UMIMatrix <- data1()
      matrixSeuratObject <- CreateSeuratObject(counts = UMIMatrix, project = "Seurat object matrix", min.cells = 0, min.features = 0)
    }
    if (VariableFeatures(matrixSeuratObject) %>% length() == 0) {
      matrixSeuratObject <- FindVariableFeatures(matrixSeuratObject, selection.method = "vst", nfeatures = 2000)
    } else {
      message("Using variable genes in object")
    }

    metadataCol <- data2()[[input$metadataCellType]]
    # use for classification of cell types
    messages <<- capture.output(
      res <- clustify(
        input = object_data(matrixSeuratObject, "data"),
        metadata = metadataCol,
        ref_mat = benchmarkRef,
        query_genes = VariableFeatures(matrixSeuratObject),
        verbose = TRUE
      ),
      type = "message"
    )
    rv$clustifym <<- messages
    
    w4$hide()
    res
  })

  output$reference <- DT::renderDataTable({
    if (rv$res_visited == 1) {
      return(df1 <- data.frame(`nodata` = rep("", 6)))
    }
    reference_matrix <- data_avg()
    rownames_to_column(as.data.frame(reference_matrix), input$metadataCellType)
  })
  
  output$clustify <- DT::renderDataTable({
    if (rv$res_visited == 1) {
      return(df1 <- data.frame(`nodata` = rep("", 6)))
    }
    res <- dataClustify()
    rownames_to_column(as.data.frame(res), input$metadataCellType)
  })
  
  corToCall <- reactive({
    res <- dataClustify()
    cor_to_call(cor_mat = res, cluster_col = input$metadataCellType)
  })
  
  output$corToCall <- DT::renderDataTable({
    corToCall()
  })
  
  # Make plots such as heat maps to compare benchmarking with clustify with actual cell types
  
  output$hmap <- renderPlot({
    if (input$metadataCellType == "") {
      return(NULL)
    }
    tmp_mat <- dataClustify()
    # could expose as an option
    cutoff_to_display <- 0.5

    if (!is.null(tmp_mat)) {
      w5$show()
    }
    tmp_mat <- tmp_mat[, colSums(tmp_mat > 0.5) > 1]
    plot_hmap(tmp_mat)
  })
  
  referenceDownload <- reactive({
    avgMatrix <- data_avg()
  })
  
  clustifyDownload <- reactive({
    clustifyMatrix <- dataClustify()
  })
  
  output$downloadReference <- downloadHandler(
    filename = function() {
      paste("reference-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(referenceDownload(), file, quote = FALSE)
    }
  )
  output$downloadClustify <- downloadHandler(
    filename = function() {
      paste("clustify-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(list(corToCall(), clustifyDownload()), file, quote = FALSE, rowNames = TRUE)
    }
  )
  
  # load example data
  observeEvent(
    input$example,
    {
      message("loading prepackaged data")
      rv$matrixloc <- list(datapath = "data/example-input/matrix.csv.gz")
      rv$metaloc <- list(datapath = "data/example-input/meta-data.csv.gz")
      updateTabItems(session, "tabs", "metadataLoad")
    }
  )

  output$clustifym <- renderUI(
    HTML(paste0(c(rv$clustifym, ""), collapse = "<br/><br/>"))
  )

  # modal for GEO id
  observeEvent(
    input$geo1 | input$geo2,
    showModal(modalDialog(
      div(id = "modalgeo",
          textInput("geoid", "query GEO id", value = rv$lastgeo),
          actionButton("geogo", "Fetch file info", icon = icon("eye"))
      ),
      easyClose = TRUE,
      fade = FALSE,
      footer = NULL
    )),
    ignoreInit = T
  )

  observeEvent(
    input$geogo,
    {
      w6$show()
      rv$lastgeo <- input$geoid
      rv$links <- list_geo(rv$lastgeo)
      message(rv$links)
      if (rv$links != "error_get") {
        rv$links2 <- rv$links %>% mutate(size = map(link, get_file_size)) %>% select(-link)
        links2 <- cbind(rv$links2,
                        button = sapply(1:nrow(rv$links), make_button("tbl1")),
                        stringsAsFactors = FALSE) %>%
          data.table::data.table()
        links2 <- links2 %>%
          DT::datatable(options = list(
            dom = "ftp",
            searchHighlight = TRUE,
            paging = TRUE,
            pageLength = 5,
            scrollY = FALSE),
            escape = ncol(links2) - 1, fillContainer = TRUE)
        
      } else {
        links2 <- data.frame(rv$links)
      }
      
      url <- str_c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", input$geoid)
      w6$hide()
      showModal(modalDialog(
        size = "l",
        div(id = "modalfiles",
            DT::renderDataTable(links2)
        ),
        easyClose = TRUE,
        fade = FALSE,
        footer = tagList(
          actionButton("geopage", label = "Go to GEO page",
                       onclick = paste0('window.open("',
                                        url,
                                        '", "_blank")'),
                       icon = icon("link")),
          actionButton("email", label = "Email author for missing data",
                       onclick = paste0('location.href="',
                                        prep_email(rv$lastgeo),
                                        '"'),
                       icon = icon("envelope-open-text")),
          actionButton("sheet", label = "Spot check for someta", 
                       icon = icon("feather-alt"))
        )
      ))
    }
  )

  observeEvent(input[["button"]], {
    w7$show()
    splitID <- strsplit(input[["button"]], "_")[[1]]
    tbl <- splitID[2]
    row <- splitID[3]
    rv$loadinglink <<- rv$links$link[as.numeric(row)]

    # if tar, read a file list
    if (str_detect(rv$links$link[as.numeric(row)], "/GSE[0-9]+_RAW.tar")) {
      fullb <- F
      previewdata <- preview_link(get_tar(rv$links$link[as.numeric(row)]))
    } else {
      fullb <- T
      previewdata <- preview_link(rv$links$link[as.numeric(row)])
    }

    if (is.null(previewdata)) {
      fullb <- F
      previewdata <- data.frame(unreadable = rep("", 4))
    } else {
      cols <- ncol(previewdata)
      previewdata <- previewdata[, 1:min(cols, 100)]
    }

    if (input[["activeTab"]] == "someta") {
      fullb <- F
    }
    
    w7$hide()
    showModal(modalDialog(
      size = "l",
      div(id = "modalback",
          title = "preview",
          DT::renderDataTable(previewdata),
          if (fullb) {
            actionButton("full", "Start full loading", icon = icon("running"))
          } else {
            disabled(actionButton("full", "Start full loading", icon = icon("running")))
          },
          actionButton("back", "Back to file list", icon = icon("step-backward"))
      ),
      easyClose = TRUE,
      fade = FALSE,
      footer = NULL
    ))
  })

  observeEvent(input$full, {
    message(rv$loadinglink)
    if (input[["activeTab"]] == "matrixLoad") {
      rv$matrixloc <- list(datapath = rv$loadinglink)
    } else if (input[["activeTab"]] == "metadataLoad") {
      rv$metaloc <- list(datapath = rv$loadinglink)
    }
    removeModal()
  })
  
  observeEvent(input$sheet, {
    showModal(modalDialog(
      size = "l",
      div(id = "modalsheet",
          title = "Please fill out",
          renderUI(h2(rv$lastgeo)),
          hr(),
          strong(materialSwitch(
            "issc",
            "   is single cell data",
            value = TRUE,
            status = "success",
            right = TRUE,
            inline = FALSE,
            width = NULL
          )),
          strong(materialSwitch(
            "hasmeta",
            "   has metadata",
            value = TRUE,
            status = "success",
            right = TRUE,
            inline = FALSE,
            width = NULL
          )),
          strong(materialSwitch(
            "hascellcol",
            "   has cell type column in metadata",
            value = TRUE,
            status = "success",
            right = TRUE,
            inline = FALSE,
            width = NULL
          )),
          textInput("comment", "", placeholder = "Additional comments")
      ),
      easyClose = TRUE,
      fade = FALSE,
      footer = actionButton("submit", "Submit", icon = icon("feather-alt"))
    ))
  })
  observeEvent(input$back, {
    links2 <- cbind(rv$links2,
                    button = sapply(1:nrow(rv$links), make_button("tbl1")),
                    stringsAsFactors = FALSE) %>%
      data.table::data.table()
    links2 <- links2 %>%
      DT::datatable(options = list(
        dom = "ftp",
        searchHighlight = TRUE,
        paging = TRUE,
        pageLength = 5,
        scrollY = FALSE),
        escape = ncol(links2)-1, fillContainer = TRUE)
    url <- str_c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", input$geoid)
    showModal(modalDialog(
      size = "l",
      div(id = "modalfiles",
          DT::renderDataTable(links2)
      ),
      easyClose = TRUE,
      fade = FALSE,
      footer = tagList(
        actionButton("geopage", label = "Go to GEO page",
                     onclick = paste0('window.open("',
                                      url,
                                      '", "_blank")'),
                     icon = icon("link")),
        actionButton("email", label = "Email author for missing data",
                     onclick = paste0('location.href="',
                                      prep_email(rv$lastgeo),
                                      '"'),
                     icon = icon("envelope-open-text")),
        actionButton("sheet", label = "Spot check for someta", 
                     icon = icon("feather-alt"))
      )
    ))
  })
  
  # upload to google sheet
  observeEvent(input$submit, {
    sheet_append(sheetid, data.frame(id = rv$lastgeo,
                                     issc = input$issc,
                                     hasmeta = input$hasmeta,
                                     hascellcol = input$hascellcol,
                                     comment = input$comment))
    # print(read_sheet(sheetid, 1))
    
    showModal(modalDialog(
      div(id = "modaldone",
          h2("Results uploaded, thank you!")
      ),
      easyClose = TRUE,
      fade = FALSE,
      footer = NULL
    ))
  })

  # disable menu at load
  addCssClass(selector = "a[data-value='clustifyres']", class = "inactiveLink")
  addCssClass(selector = "ul li:eq(4)", class = "inactiveLink")
  
  addCssClass(selector = "a[data-value='blank']", class = "inactiveLink")
  addCssClass(selector = "ul li:eq(5)", class = "inactiveLink")

  # check if data is loaded
  observeEvent((!is.null(data1())) + (!is.null(data2())) + (!is.null(data3())) +
                 (!is.null(input$metadataCellType)) +
                 (input$metadataCellType != ""), {
                   if ((!is.null(data1())) + (!is.null(data2())) + (!is.null(data3())) + 
                       (!is.null(input$metadataCellType)) + 
                       (input$metadataCellType != "") == 5) {
                     removeCssClass(selector = "a[data-value='clustifyres']", class = "inactiveLink")
                     removeClass(selector = "ul li:eq(4)", class = "inactiveLink")
                   }
                 })
  observeEvent(data1(), {
    if (!is.null(data1())) {
      addCssClass(selector = "a[data-value='matrixLoad']", class = "doneLink")
      addClass(selector = "ul li:eq(1)", class = "doneLink")
    }
  })

  observeEvent(input$metadataCellType, {
    if (input$metadataCellType != "") {
      addCssClass(selector = "a[data-value='metadataLoad']", class = "doneLink")
      addClass(selector = "ul li:eq(2)", class = "doneLink")
    }
  })

  observeEvent(input[["activeTab"]], {
    if (input[["activeTab"]] == "clusterRef") {
      rv$ref_visited <<- 1
    } else if (input[["activeTab"]] == "clustifyres") {
      if (rv$res_visited == 0)
        rv$res_visited <<- 1
    }
  })

  observeEvent(rv$ref_visited, {
    if (rv$ref_visited == 1 & !is.null(data3())) {
      addCssClass(selector = "a[data-value='clusterRef']", class = "doneLink")
      addClass(selector = "ul li:eq(3)", class = "doneLink")
    }
  })

  observeEvent(rv$res_visited, {
    if (rv$res_visited == 1) {
      rv$res_visited <- 2
    }
  }, ignoreInit = FALSE)

  observeEvent(rv$ref_link, {
    runjs(paste0("document.getElementById('ref_linkgo').onclick = function() {
           window.open('", rv$ref_link, "', '_blank');};"))
  })
  
  output$someta <- DT::renderDataTable({
    as.data.table(someta %>% select(-geo, -pubmed, -pubmed_id), rownames = FALSE)
  }, filter = "top", selection=list(mode="single", target="row"),
  rownames = FALSE, options = list(autoWidth = TRUE,
                                                      columnDefs = list(
    list(width = '200px', targets = c(0:6)), list(
    targets = c(3, 4,5),
    render = JS(
      "function(data, type, row, meta) {",
      "return type === 'display' && data.length > 100 ?",
      "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
      "}")
  ))))
  
  observeEvent(input$someta_cell_clicked, {
    if (length(input$someta_cell_clicked) != 0) {
      sel <- input$someta_cell_clicked
      rv$lastgeo <- someta$id[sel$row]
      w6$show()
      rv$links <- list_geo(rv$lastgeo)
      message(rv$links)
      if (rv$links != "error_get") {
        rv$links2 <- rv$links %>% mutate(size = map(link, get_file_size)) %>% select(-link)
        links2 <- cbind(rv$links2,
                        button = sapply(1:nrow(rv$links), make_button("tbl1")),
                        stringsAsFactors = FALSE) %>%
          data.table::data.table()
        links2 <- links2 %>%
          DT::datatable(options = list(
            dom = "ftp",
            searchHighlight = TRUE,
            paging = TRUE,
            pageLength = 5,
            scrollY = FALSE),
            escape = ncol(links2) - 1, fillContainer = TRUE)
        
      } else {
        links2 <- data.frame(rv$links)
      }
      
      url <- str_c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", input$geoid)
      w6$hide()
      showModal(modalDialog(
        size = "l",
        div(id = "modalfiles",
            DT::renderDataTable(links2)
        ),
        easyClose = TRUE,
        fade = FALSE,
        footer = tagList(
          actionButton("geopage", label = "Go to GEO page",
                       onclick = paste0('window.open("',
                                        url,
                                        '", "_blank")'),
                       icon = icon("link")),
          actionButton("email", label = "Email author for missing data",
                       onclick = paste0('location.href="',
                                        prep_email(rv$lastgeo),
                                        '"'),
                       icon = icon("envelope-open-text")),
          actionButton("sheet", label = "Spot check for someta", 
                       icon = icon("feather-alt"))
        )
      ))
    }
  })
}