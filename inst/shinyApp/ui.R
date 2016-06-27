# App ui

headerbar <- dashboardHeader(
  title = "ALTRE Analysis Workflow",
  titleWidth = 270,
  dropdownMenu(
    type = "notifications",
    notificationItem(
      text = "Plots take some time to display",
      icon("truck"),
      status = "warning"
      )
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
      "Categorize Altered Regions",
      icon = icon("bullseye"),
      tabName = "cataltered",
      badgeLabel = "step 6",
      badgeColor = "green"
    ),
    menuItem(
      "Compare Altered Regions",
      icon = icon("balance-scale"),
      tabName = "compare",
      badgeLabel = "step 7",
      badgeColor = "green"
    ),
    menuItem(
      "Find Enriched Pathways",
      icon = icon("gears"),
      tabName = "pathways",
      badgeLabel = "step 8",
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
                       width = 12,
                       includeMarkdown("include.md")
                       )
                     )
            ),
    tabItem(tabName = "loaddata",
            fluidRow(
              HTML("<div class='col-sm-6' style='min-width:
                   700px !important;'>"),
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
              HTML("</div>"),
              infoBoxOutput("statusbox1", width = 6)
              )
            ),
    tabItem(tabName = "definerep",
            fluidRow(
              HTML("<div class='col-sm-5' style='min-width:
                   600px !important;'>"),
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
              infoBoxOutput("statusbox2", width = 7),
              HTML("<div class='col-sm-7' style='min-width:
                   500px !important;'>"),
              box(
                title = "Barplot",
                width = NULL,
                solidHeader = TRUE,
                plotOutput('barplot')
              ),
              HTML("</div>")
              )
            ),
    tabItem(tabName = "combine",
            fluidRow(
              HTML("<div class='col-sm-4' style='min-width:
                   300px !important;'>"),
              box(
                width = NULL,
                solidHeader = TRUE,
                title = "Combine and Annotate Peaks",
                actionButton("buttonannot", strong("Combine and Annotate")),
                radioButtons(
                  "mergeradio",
                  label = h4("Merge"),
                  choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                  selected = "TRUE"
                ),
                sliderInput(
                  "distTSS",
                  label = h4("Distance from TSS"),
                  min = 0,
                  max = 3000,
                  value = 1500
                ),
                hr(),
                radioButtons(
                  "regionradio",
                  label = h4("Region specific merging?"),
                  choices = list("TRUE" = "TRUE", "FALSE" = "FALSE")
                ),
                conditionalPanel("input.regionradio == 'FALSE'",
                                 sliderInput(
                                   "dist",
                                   label = h4("Merge promoters and enhancers if
                                              distance is less than "),
                                   min = 0,
                                   max = 3000,
                                   value = 0
                                   )
                                 ),
                conditionalPanel("input.regionradio == 'TRUE'",
                                 sliderInput(
                                   "distenh",
                                   label = h4("Merge enhancers distance
                                              threshold"),
                                   min = 0,
                                   max = 3000,
                                   value = 1500
                                   ),
                                 sliderInput(
                                   "distprom",
                                   label = h4("Merge promoters distance
                                              threshold"),
                                   min = 0,
                                   max = 3000,
                                   value = 1000
                                   )
                                 )
                ),
              HTML("</div>"),
              infoBoxOutput("statusbox3", width = 7),
              HTML("<div class='col-sm-7' style='min-width:
                   500px !important;'>"),
              box(
                title = "Barplot",
                width = NULL,
                solidHeader = TRUE,
                plotOutput('annotatebarplot'),
                dataTableOutput("table3")
              ),
              HTML("</div>")
            )),
    tabItem(tabName = "retrieve",
            fluidRow(
              HTML("<div class='col-sm-4' style='min-width:
                   300px !important;'>"),
              box(
                width = NULL,
                title = "Retrieve Read Counts in Annotated Regions",
                actionButton("buttoncounts", strong("Retrieve Counts")),
                hr(),
                uiOutput("chooseref")
              ),
              HTML("</div>"),
              infoBoxOutput("statusbox4", width = 7),
              HTML("<div class='col-sm-7' style='min-width:
                   500px !important;'>"),
              box(
                width = NULL,
                title = "Density Plot",
                plotOutput('densityplot')
              ),
              HTML("</div>")
              )
            ),
    tabItem(tabName = "definealtered",
            fluidRow(
              HTML("<div class='col-sm-4' style='min-width:
                   300px !important;'>"),
              box(
                width = NULL,
                title = "Define Altered Regions",
                actionButton("buttondefine", strong("Define Altered Regions")),
                hr(),
                sliderInput(
                  "alpha",
                  label = h4("pvalue cutoff Cutoff"),
                  min = 0,
                  max = 1,
                  value = 0.01
                ),
                sliderInput(
                  "lfcThreshold",
                  label = h4("log2fold change cutoff"),
                  min = 0,
                  max = 5,
                  value = 1,
                  step = 0.1
                )
              ),
              HTML("</div>"),
              infoBoxOutput("statusbox5", width = 7),
              HTML("<div class='col-sm-7' style='min-width:
                   500px !important;'>"),
              HTML("</div>")
              )
            ),
    tabItem(tabName = "cataltered",
            fluidRow(
              HTML("<div class='col-sm-4' style='min-width:
                   300px !important;'>"),
              box(
                width = NULL,
                title = "Categorize Altered Regions",
                actionButton("buttoncat", strong("Categorize Altered Regions")),
                hr(),
                sliderInput(
                  "lfcSpecific",
                  label = h4("log2fold change cutoff for specific
                             enhancers/promoters"),
                  min = 0,
                  max = 5,
                  value = 1.5,
                  step=0.1
                ),
                sliderInput(
                  "lfcShared",
                  label = h4("log2fold change cutoff for shared
                             enhancers/promoters"),
                  min = 0,
                  max = 5,
                  value = 1.2,
                  step = 0.1
                ),
                sliderInput(
                  "pvalueSpecific",
                  label = h4("pvalue cutoff for specific
                             enhancers/promoters"),
                  min = 0,
                  max = 1,
                  value = 0.01
                ),
                sliderInput(
                  "pvalueShared",
                  label = h4("pvalue cutoff for shared enhancers/promoters"),
                  min = 0,
                  max = 1,
                  value = 0.05
                )
              ),
              HTML("</div>"),
              infoBoxOutput("statusbox5b", width=7),
              HTML("<div class='col-sm-7' style='min-width:
                   500px !important;'>"),
              box(
                width = NULL,
                title = "Volcano Plot",
                plotOutput('volcanoplot'),
                dataTableOutput("table4")
              ),
              HTML("</div>")
              )
  ),
    tabItem(tabName = "compare",
            fluidRow(
              HTML("<div class='col-sm-4' style='min-width:
                   300px !important;'>"),
              box(
                width = NULL,
                title = "Compare methods of identifying altered
                regulatory regions",
                actionButton("buttoncompare", strong("Compare Methods")),
                hr(),
                dataTableOutput("table5")
              ),
              HTML("</div>"),
              infoBoxOutput("statusbox8", width = 7),
              HTML("<div class='col-sm-7' style='min-width:
                   500px !important;'>"),
              box(
                width = NULL,
                title = "Venn Plot",
                plotOutput('vennplot')
              ),
              HTML("</div>")
              )
            ),
    tabItem(tabName = "pathways",
            fluidRow(
              infoBoxOutput("statusbox6", width = 6),
              infoBoxOutput("statusbox7", width = 6)
            ),
            fluidRow(
              box(
                title = "Pathway enrichment MF",
                actionButton("buttonpathwayMF",
                             strong("Run Pathway Enrichment MF")),
                hr(),
                sliderInput(
                  "pathpvaluecutoffMF",
                  label = h4("pvalue cutoff"),
                  min = 0,
                  max = 1,
                  value = 0.01
                  ),
                hr(),
                plotOutput('heatplotMF')
                ),
              box(
                title = "Pathway enrichment BP",
                actionButton("buttonpathwayBP",
                             strong("Run Pathway Enrichment BP")),
                hr(),
                sliderInput(
                  "pathpvaluecutoffBP",
                  label = h4("pvalue cutoff"),
                  min = 0,
                  max = 1,
                  value = 0.01
                  ),
                hr(),
                plotOutput('heatplotBP')
                )
              )
            )
    )
)

shinyUI(
  dashboardPage(
    headerbar,
    sidebar,
    body
    )
  )
