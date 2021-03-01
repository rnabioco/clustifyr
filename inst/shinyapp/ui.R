# Define UI for data upload app ----
ui <- dashboardPage(
  title = "clustifyr app",
  skin = "green",
  dashboardHeader(title = div(tags$a(href='https://github.com/rnabioco/clustifyr',
                                 tags$img(src='logo.png', width="18%")),
                                 "clustifyr Shiny app")),
  dashboardSidebar(
    # tags$img(
    #   src = 'logo.png',
    #   width = "1000%",
    #   style = 'position: fixed; bottom: 0;right: 0;'
    # ),
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
    shinyDashboardThemes(theme = "poor_mans_flatly"),
    tabItems(
      tabItem(
        tabName = "dashboard",
        # js stuff ----
        useShinyjs(),
        use_bs_tooltip(),
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
                color: green;
            }'
        ))),
        tags$head(tags$style(HTML('
            .skin-green .sidebar .doneLink.active > a {
                color: green;
                border-left-color: green;
            }'
        ))),
        tags$head(tags$style(HTML('
            .skin-green .sidebar .doneLink:hover {
                color: green;
                border-left-color: green;
            }'
        ))),

        # waiter stuff ----
        use_waiter(),

        # load example data ----
        actionButton("example",
                     "Load example data",
                     icon = icon("space-shuttle")
        ) %>%
          bs_embed_tooltip("Use example data from GSE113049 for walkthrough", placement = "right"),
        
        # Horizontal line ----
        tags$hr(),
        
        # readme
        includeMarkdown("README.md"),
        
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
        h2("Load Counts (raw or normalized) Matrix"),
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
                    ".rdata",
                    ".gz"
                  )
        ) %>% 
          bs_embed_tooltip("Accepted file types: plain or gz text files, and rds/rda/rdata for Seurat/SCE",
                           placement = "right"),

        # GEO id load ----
        actionButton("geo1",
                     "or from GEO id",
                     icon = icon("search")) %>% 
          bs_embed_tooltip("Alternatively, load data directly from GEO with ID", placement = "right"),

        actionButton("matrixPopup", "Display UMI Matrix in popup"),
        tags$hr(),
        DTOutput("contents1"), # UMI Count Matrix
        tags$hr()
      ),
      tabItem(
        tabName = "metadataLoad",
        h2("Load Metadata"),
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
                    ".rdata",
                    ".gz"
                  )
        ) %>% 
          bs_embed_tooltip("Accepted file types: plain or gz text files, and rds/rda/rdata for Seurat/SCE",
                           placement = "right"),

        # GEO id load ----
        actionButton("geo2",
                     "or from GEO id",
                     icon = icon("search")) %>% 
          bs_embed_tooltip("Alternatively, load data directly from GEO with ID", placement = "right"),

        actionButton("metadataPopup", "Display Metadata table in popup"),
        h2("Choose column in metadata with cluster info"),
        selectInput("metadataCellType", NULL,
                    choice = list("")
        ) %>% 
          bs_embed_tooltip("Select from dropdown, or click on preview column below", placement = "right"),
        tags$hr(),
        DTOutput('contents2'),
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
        ) %>% 
          bs_embed_tooltip("Select from pre-built references in clustifyrdatahub", placement = "right"),
        actionButton("ref_linkgo",
                     label = "Go to original source",
                     icon = icon("link")) %>% 
          bs_embed_tooltip("For more info on the reference datasets", placement = "right"),
        tags$hr(),
        h2("Or load reference table"),
        fileInput("file3", "Choose Reference Average Expression File",
                  multiple = FALSE,
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    ".tsv",
                    ".xlsx",
                    ".gz"
                  )
        ) %>% 
          bs_embed_tooltip("Alternatively, choose local file of average gene expression by cell type", placement = "right"),
        tags$hr(),
        uiOutput("ref_summary"),
        DTOutput("contents3", height = "300px"),
        tags$hr()
      ),
      tabItem(
        tabName = "clustifyres",
        box(id = "box_clustifym",
            collapsible = TRUE,
            collapsed = TRUE,
            solidHeader = TRUE,
            status = "info",
            title = "clustifyr messages",
            htmlOutput("clustifym")) %>% 
          bs_embed_tooltip("Console messages from run", "right"),
        div(style="display:inline-block",
          downloadButton("downloadReference", "Download average expression matrix") %>% 
            bs_embed_tooltip("Note: this is the same file needed to use the current dataset as reference in the future", "bottom"),
          downloadButton("downloadClustify", "Download clustify results matrices") %>% 
            bs_embed_tooltip("Cell type inference as sheet1, correlation matrix as sheet2, of xlsx", "bottom")
        ),
        tags$hr(),
        # actionButton("uploadClustify", "Upload reference matrix"),
        h2("Average Expression Matrix"),
        DT::dataTableOutput("reference", height = "300px"), # Reference Matrix
        tags$hr(),
        h2("Ranked Correlation Matrix"), 
        DT::dataTableOutput("clustify", height = "300px"), # Clustify Matrix
        tags$hr(),
        h2("Cell Type Inference Results"),
        DT::dataTableOutput("corToCall", height = "300px"),
        tags$hr(),
        plotOutput("hmap", height = "900px"),
        tags$hr(),
      ),
      tabItem(
        tabName = "someta",
        h2("GEO scRNA-seq records, hover to see truncated text, click to preview files"),
        DT::dataTableOutput("someta", height = "1000px")
      )
    )
  )
)