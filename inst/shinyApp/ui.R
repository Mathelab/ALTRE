# App ui

headerbar <- dashboardHeader(
  title = "ALTRE Analysis Workflow",
  titleWidth = 270,
  dropdownMenu(
    type = "notifications",
    notificationItem(text = "Plots take some time to display",
                     icon("truck"),
                     status = "warning")
  )
)

sidebar <- dashboardSidebar(
  width = 270,
  sidebarMenu(
    menuItem("About",
             tabName = "about",
             icon = icon("info")),
    menuItem(
      "Load Data",
      tabName = "loaddata",
      icon = icon("folder-open"),
      badgeLabel = "step 1"
    ),
    menuItem(
      "Define Consensus Peaks",
      tabName = "definerep",
      icon = icon("bullseye"),
      badgeLabel = "step 2"
    ),
    menuItem(
      "Annotate Peaks",
      tabName = "combine",
      icon = icon("picture-o"),
      badgeLabel = "step 3"
    ),
    menuItem(
      "Retrieve Read Counts",
      icon = icon("bolt"),
      tabName = "retrieve",
      badgeLabel = "step 4",
      badgeColor = "green"
    ),
    menuItem(
      "Define Altered Regions",
      icon = icon("bullseye"),
      tabName = "definealtered",
      badgeLabel = "step 5",
      badgeColor = "green"
    ),
    menuItem(
      "Find Enriched Pathways",
      icon = icon("gears"),
      tabName = "pathways",
      badgeLabel = "step 6",
      badgeColor = "green"
    ),
    menuItem(
      "Compare Altered Regions",
      icon = icon("balance-scale"),
      tabName = "compare",
      badgeLabel = "step 7",
      badgeColor = "green"
    )
  )
)

body <- dashboardBody(
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
  ),

  tabItems(
    tabItem(tabName = "about",
            tabPanel("About",
                     box(
                       width = 12, includeMarkdown("include.md")
                     ))),
    tabItem(tabName = "loaddata",
            fluidRow(
              HTML("<div class='col-sm-6' style='min-width: 700px !important;'>"),
              box(
                title = "Load Datapaths File",
                width = NULL,
                solidHeader = TRUE,
                fileInput(
                  "file",
                  accept = c('text/csv',
                             'text/comma-separated-values,text/plain',
                             '.csv'),
                  "Provide CSV File to Load Data:"
                ),
                dataTableOutput("table1")
              ),
              HTML("</div>")
            ))
    ,
    tabItem(tabName = "definerep",
            fluidRow(
              HTML("<div class='col-sm-4' style='min-width: 600px !important;'>"),
              box(
                title = "Load and Merge Peak Files",
                width = NULL,
                solidHeader = TRUE,
                numericInput(
                  "numOverlap",
                  "Minimum Number of Overlaps to Determine Consensus Region ",
                  2,
                  min = 2,
                  max = 10
                ),
                actionButton("buttonmerge", strong("Merge Replicates")),
                hr(),
                dataTableOutput("table2")
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-6' style='min-width: 500px !important;'>"),
              box(
                title = "Barplot",
                width = NULL,
                solidHeader = TRUE,
                plotOutput('barplot')
              ),
              HTML("</div>")
            )),
    tabItem(tabName = "combine",
            fluidRow(
              box(
                title = "Combine and Annotate Peaks",
                actionButton("buttonannot", strong("Combine and Annotate"))
              ),
              box(
                title = "Table",
                dataTableOutput("table3")
            )
            )
            ),
    tabItem(tabName = "retrieve",
            fluidRow(
              box(title = "Retrieve Read Counts in Annotated Regions",
                  actionButton(
                    "buttoncounts", strong("Retrieve Counts")
                  )),

              box(title = "Density Plot",
                  plotOutput('densityplot'))
            )),
    tabItem(tabName = "definealtered",
            fluidRow(
              box(title = "Define Altered Regions",
                  actionButton("buttondefine", strong("Define"))),
              box(title = "Volcano Plot",
                  plotOutput('volcanoplot'))
            )),
    tabItem(tabName = "pathways",
            fluidRow(
              box(title = "Pathway enrichment",
                  actionButton(
                    "buttonpathway", strong("Run Pathway Enrichment")
                  )),
              box(title = "Enrichment Heatmap",
                  plotOutput('heatplot'))
            ))

  )
)

shinyUI(dashboardPage(headerbar,
                      sidebar,
                      body))
