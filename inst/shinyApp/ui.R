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
                title = "Load File",
                width = NULL,
                solidHeader = TRUE,
                fileInput(
                  "file",
                  accept = c('text/csv',
                             'text/comma-separated-values,text/plain',
                             '.csv'),
                  "Locate CSV Metadata File"
                  ),
                dataTableOutput("table1")
                ),
              HTML("</div>"),
              infoBoxOutput("statusbox1", width = 6)
              )
            ),
    tabItem(tabName = "definerep",
            fluidRow(
              HTML("<div class='col-sm-4' style='min-width:
                   550px !important;'>"),
              box(
                title = strong("Load and Merge Annotation Files") ,
                width = NULL,
                solidHeader = TRUE,
                h5("This step does the following: "),
                tags$ol(
                  tags$li("Loads biosample annotation files of DNase I hypersensitive
                   sites (i.e peaks)"),
                  tags$li(" Merges the bioreplicates of each biosample in order
                          to determine the consensus peaks, defined as genomic
                          regions that overlap in at least N of the bioreplicates.")
                  ),
                hr(),
                numericInput(
                  "numOverlap",
                  "Choose the Minimum Number of Overlaping Genomic Regions N",
                  2,
                  min = 2,
                  max = 10
                ),
                actionButton("buttonmerge", strong("Load files then Merge Replicates")),
                hr(),
                dataTableOutput("table2")
              ),
              HTML("</div>"),
              infoBoxOutput("statusbox2", width = 7),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
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
                h5("Combine and annotate peaks from different sample types.
                   Optionally merge nearby regions."),
                actionButton("buttonannot", strong("Combine and Annotate")),
                radioButtons(
                  "mergeradio",
                  label = h5("Merge"),
                  choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                  selected = "TRUE"
                ),
                sliderInput(
                  "distTSS",
                  label = h5("Distance from TSS"),
                  min = 0,
                  max = 3000,
                  value = 1500
                ),
                hr(),
                radioButtons(
                  "regionradio",
                  label = h5("Region specific merging?"),
                  choices = list("TRUE" = "TRUE", "FALSE" = "FALSE")
                ),
                conditionalPanel("input.regionradio == 'FALSE'",
                                 sliderInput(
                                   "dist",
                                   label = h5("Merge promoters and enhancers if
                                              distance is less than "),
                                   min = 0,
                                   max = 3000,
                                   value = 0
                                   )
                                 ),
                conditionalPanel("input.regionradio == 'TRUE'",
                                 sliderInput(
                                   "distenh",
                                   label = h5("Merge enhancers distance
                                              threshold"),
                                   min = 0,
                                   max = 3000,
                                   value = 1500
                                   ),
                                 sliderInput(
                                   "distprom",
                                   label = h5("Merge promoters distance
                                              threshold"),
                                   min = 0,
                                   max = 3000,
                                   value = 1000
                                   )
                                 )
                ),
              HTML("</div>"),
              infoBoxOutput("statusbox3", width = 8),
              HTML("<div class='col-sm-8' style='min-width:
                   500px !important;'>"),
              box(
                title = "Barplot",
                width = NULL,
                solidHeader = TRUE,
                plotOutput('annotatebarplot'),
                hr(),
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
                title = "Retrieve Read Counts",
                h5("Counts the number of reads in each regulatory region
                   for each sample type."),
                actionButton("buttoncounts", strong("Retrieve Counts")),
                hr(),
                uiOutput("chooseref"),
                hr(),
                radioButtons(
                  "chromradio",
                  label = h5("Restrict Analysis to a Single Chromosome?"),
                  choices = list("FALSE" = "FALSE","TRUE" = "TRUE"),
                  selected = "FALSE"
                ),
                conditionalPanel("input.chromradio == 'TRUE'",
                uiOutput("chooseChrom")
                )
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
                h5(" Determines which regulatory regions are signifigantly
                   altered between sample types."),
                actionButton("buttondefine", strong("Define Altered Regions")),
                hr(),
                sliderInput(
                  "alpha",
                  label = h5("pvalue cutoff"),
                  min = 0,
                  max = 1,
                  value = 0.01
                ),
                sliderInput(
                  "lfcThreshold",
                  label = h5("log2fold change cutoff"),
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
                title = "Categorize Regions",
                h5("Categorize altered regulatory regions as experiment-specific,
                reference-specific, or shared."),
                actionButton("buttoncat", strong("Categorize Altered Regions")),
                hr(),
                h4("Parameters that define cell-type specific regulatory regions."),
                sliderInput(
                  "lfcSpecific",
                  label = h5("log2fold change cutoff for specific
                             enhancers/promoters"),
                  min = 0,
                  max = 5,
                  value = 1.5,
                  step = 0.1
                ),
                sliderInput(
                  "pvalueSpecific",
                  label = h5("pvalue cutoff for specific
                             enhancers/promoters"),
                  min = 0,
                  max = 1,
                  value = 0.01
                ),
                hr(),
                h4("Parameters that define shared regulatory regions."),
                sliderInput(
                  "lfcShared",
                  label = h5("log2fold change cutoff for shared
                             enhancers/promoters"),
                  min = 0,
                  max = 5,
                  value = 1.2,
                  step = 0.1
                ),
                sliderInput(
                  "pvalueShared",
                  label = h5("pvalue cutoff for shared enhancers/promoters"),
                  min = 0,
                  max = 1,
                  value = 0.05
                )
                #,
                # hr(),
                # conditionalPanel("input.buttoncat > 0",
                #                  downloadButton("downloadData",
                #                               strong("Download Track")
                #                               )
                #                  )
              ),
              HTML("</div>"),
              infoBoxOutput("statusbox5b", width = 7),
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
                title = "Compare Methods",
                h5("Compare two methods of identifying altered regulatory regions,
                   one based on peak intensity, the other on peak presence."),
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
    tabItem(
      tabName ="pathways",
      tabBox(
        title = "Pathway Enrichment Analysis",
        width = 12,
        tabPanel("Pathway Enrichment for Molecular Function",
                 fluidRow(
                   infoBoxOutput("statusbox6", width = 11)
                   ),
                 fluidRow(
                   box(
                     title = "Pathway Enrichment for Molecular Function",
                     width = 11,
                     h5("Determine which pathways are overrepresented in altered
                        promoters and enhancers."),
                     actionButton("buttonpathwayMF",
                                  strong("Run MF Pathway Enrichment")),
                     hr(),
                     sliderInput(
                       "pathpvaluecutoffMF",
                       label = h5("pvalue cutoff"),
                       min = 0,
                       max = 1,
                       value = 0.01
                     ),
                     hr(),
                     plotOutput('heatplotMF')
                     )
                   )
                 ),
        tabPanel("Pathway Enrichment for Biological Process",
               fluidRow(
                 infoBoxOutput("statusbox7", width = 11)
                 ),
               fluidRow(
                 box(
                   title = "Pathway Enrichment for Biological Process",
                   width = 11,
                   h5("Determine which pathways are overrepresented in altered promoters
                      and enhancers."),
                   actionButton("buttonpathwayBP",
                                strong("Run BP Pathway Enrichment")),
                   hr(),
                   sliderInput(
                     "pathpvaluecutoffBP",
                     label = h5("pvalue cutoff"),
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
  )
)

shinyUI(
  dashboardPage(
    headerbar,
    sidebar,
    body
    )
  )
