# Define UI for data upload app ----
ui <- dashboardPage(
  skin = "green",
  dashboardHeader(title = "clustifyr RShiny app"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Load Matrix", tabName = "matrixLoad", icon = icon("th")),
      menuItem("Load Metadata", tabName = "metadataLoad", icon = icon("list")),
      menuItem("Load Reference", tabName = "clusterRef", icon = icon("database")),
      menuItem("Cell type inference", tabName = "clustifyres", icon = icon("calculator")),
      menuItem("", tabName = "blank"),
      menuItem("Explore GEO", tabName = "someta", icon = icon("hdd"))
    )
  ),
  dashboardBody(
    #shinyDashboardThemes(theme = "flat_red"),
    tabItems(
      tabItem(
        tabName = "dashboard",
        # js stuff ----
        useShinyjs(),
        tags$head(
          tags$script(HTML(js2))
        ),
        tags$head(tags$style(HTML('
            .skin-green .sidebar .inactiveLink {
                color: black;
                opacity : 25%;
            }'
        ))),
        tags$head(tags$style(".inactiveLink {
                           pointer-events: none;
                           cursor: not-allowed;
                           }")),
        tags$head(tags$style(HTML('
            .skin-green .sidebar .doneLink {
                color: blue;
            }'
        ))),
        tags$head(tags$style(HTML('
            .skin-green .sidebar .doneLink.active > a {
                color: blue;
                border-left-color: blue;
            }'
        ))),
        tags$head(tags$style(HTML('
            .skin-green .sidebar .doneLink:hover {
                color: blue;
                border-left-color: blue;
            }'
        ))),

        # waiter stuff ----
        use_waiter(),

        # load example data ----
        actionButton("example",
                     "load example data",
                     icon = icon("space-shuttle")
        ),
        
        # Input: Checkbox if file has header ----
        checkboxInput("header", "Header", TRUE),

        # Horizontal line ----
        tags$hr(),

        # Input: Select separator ----
        radioButtons("sepMat", "Separator - Matrix",
                     choices = c(
                       Comma = ",",
                       Semicolon = ";",
                       Tab = "\t"
                     ),
                     selected = ","
        ),

        radioButtons("sepMeta", "Separator - Metadata",
                     choices = c(
                       Comma = ",",
                       Semicolon = ";",
                       Tab = "\t"
                     ),
                     selected = ","
        ),

        # Input: Select number of rows to display ----
        radioButtons("dispMat", "Display - Matrix",
                     choices = c(
                       Head = "head",
                       All = "all"
                     ),
                     selected = "head"
        ),
        radioButtons("dispMeta", "Display - Metadata",
                     choices = c(
                       Head = "head",
                       All = "all"
                     ),
                     selected = "head"
        )
      ),
      tabItem(
        tabName = "matrixLoad",
        h2("Load UMI Counts Matrix"),
        # Input: Select a file ----
        fileInput("file1", "Choose Matrix File",
                  multiple = TRUE,
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    ".xlsx",
                    ".tsv",
                    ".rds",
                    ".rda",
                    ".rdata"
                  )
        ),

        # GEO id load ----
        actionButton("geo1",
                     "from GEO id",
                     icon = icon("search")),

        actionButton("matrixPopup", "Display UMI Matrix in popup"),
        DTOutput("contents1"), # UMI Count Matrix
        tags$hr()
      ),
      tabItem(
        tabName = "metadataLoad",
        h2("Load Metadata table"),
        fileInput("file2", "Choose Metadata File",
                  multiple = FALSE,
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    ".xlsx",
                    ".tsv",
                    ".rds",
                    ".rda",
                    ".rdata"
                  )
        ),

        # GEO id load ----
        actionButton("geo2",
                     "from GEO id",
                     icon = icon("search")),

        actionButton("metadataPopup", "Display Metadata table in popup"),
        h2("Choose column in metadata with cluster information"),
        selectInput("metadataCellType", "Cell Type Metadata Column:",
                    choice = list("")
        ),
        fluidRow(column(12, DTOutput('contents2'))),
        #DT::dataTableOutput("contents2"), # Metadata table
        tags$hr(),
        uiOutput("colclicked")
      ),
      tabItem(
        tabName = "clusterRef",
        h2("Choose built-in reference dataset or upload average expression matrix"),
        selectInput("dataHubReference", "ClustifyrDataHub Reference:",
                    choices = list(
                      "ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs",
                      "ref_mouse.rnaseq", "ref_moca_main", "ref_immgen", "ref_hema_microarray",
                      "ref_cortex_dev", "ref_pan_indrop", "ref_pan_smartseq2",
                      "ref_mouse_atlas"
                    )
        ),
        actionButton("ref_linkgo",
                     label = "Go to original source",
                     icon = icon("link")),
        tags$hr(),
        h2("Or load reference table"),
        fileInput("file3", "Choose Reference Average Expression File",
                  multiple = FALSE,
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    ".xlsx",
                    ".tsv",
                    ".rds",
                    ".rda"
                  )
        ),
        DTOutput("contents3", height = "300px"),
        uiOutput("ref_summary")
      ),
      tabItem(
        tabName = "clustifyres",
        box(id = "box_clustifym",
            collapsible = TRUE,
            collapsed = TRUE,
            solidHeader = TRUE,
            status = "info",
            title = "clustifyr messages",
            htmlOutput("clustifym")),
        downloadButton("downloadReference", "Download reference matrix"),
        downloadButton("downloadClustify", "Download clustify matrix"),
        # actionButton("uploadClustify", "Upload reference matrix"),
        h2("average matrix"),
        DT::dataTableOutput("reference", height = "300px"), # Reference Matrix
        tags$hr(),
        h2("ranked correlation matrix"), 
        DT::dataTableOutput("clustify", height = "300px"), # Clustify Matrix
        tags$hr(),
        h2("cell type results"),
        DT::dataTableOutput("corToCall", height = "300px"),
        tags$hr(),
        plotOutput("hmap", height = "900px")
      ),
      tabItem(
        tabName = "someta",
        DT::dataTableOutput("someta", height = "1000px")
      )
    )
  )
)