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
              HTML("<div class='col-sm-4' style='min-width:
                   600px !important;'>"),
              box(
                title = strong("Load Metadata Spreadsheet"),
                width = NULL,
                solidHeader = TRUE,
                h5("This step does the following: "),
                tags$ul(
                  tags$li("Loads a metadata spreadsheet with a csv file extention."),
                  tags$li(" Prints out the contents of entire the file except for the
                        first column, which contains the data filepaths")
                ),
                hr(),
                fileInput(
                  "file",
                  accept = c('text/csv',
                             'text/comma-separated-values,text/plain',
                             '.csv'),
                  "Lood CSV File"
                  ),
                hr(),
                dataTableOutput("table1")
                ),
              HTML("</div>"),
              infoBoxOutput("statusbox1", width = 7)
              )
            ),
    tabItem(tabName = "definerep",
            fluidRow(
              HTML("<div class='col-sm-4' style='min-width:
                   400px !important;'>"),
              box(
                title = strong("Load and Merge Annotation Files") ,
                width = NULL,
                solidHeader = TRUE,
                h5("This step does the following: "),
                tags$ul(
                  tags$li("Loads biosample annotation files of DNase I hypersensitive
                   sites (i.e. peaks)."),
                  tags$li(" Merges the bioreplicates of each biosample in order
                          to determine the consensus peaks, defined as genomic
                          regions that overlap in at least N of the bioreplicates."),
                  tags$li(" Outputs a barplot summary of the number of peaks in
                          each replicate and in their merged consensus.")
                  ),
                hr(),
                numericInput(
                  "numOverlap",
                  "Choose the Minimum Number of Overlaping Genomic Regions N",
                  2,
                  min = 2,
                  max = 10
                ),
                hr(),
                actionButton("buttonmerge", strong("Load Files then Merge Replicates")),
                hr(),
                dataTableOutput("table2")
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              infoBoxOutput("statusbox2", width = NULL),
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
                   400px !important;'>"),
              box(
                width = NULL,
                solidHeader = TRUE,
                title = strong("Combine and Annotate Peaks"),
                h5("This step does the following: "),
                tags$ul(
                  tags$li("Combines peaks from different sample types,
                   optionally merging nearby regions."),
                  tags$li(" Annotates genomic regions with type specificity based
                          on whether each region is a candidate promoter or enhancer
                          as determined by user-defined distance from a
                          transcription start site.")
                ),
                hr(),
                h4(" Select Parameters"),
                br(),
                radioButtons(
                  "mergeradio",
                  label = strong("Merge?"),
                  choices = list("TRUE" = "TRUE", "FALSE" = "FALSE"),
                  selected = "TRUE"
                ),
                hr(),
                sliderInput(
                  "distTSS",
                  label = strong("Select distance from TSS"),
                  min = 0,
                  max = 3000,
                  value = 1500
                ),
                hr(),
                radioButtons(
                  "regionradio",
                  label = strong("Region specific merging?"),
                  choices = list("TRUE" = "TRUE", "FALSE" = "FALSE")
                ),
                conditionalPanel("input.regionradio == 'FALSE'",
                                 sliderInput(
                                   "dist",
                                   label = strong("Select upper threshold distance
                                                  for merging promoters and enhancers"),
                                   min = 0,
                                   max = 3000,
                                   value = 0
                                   )
                                 ),
                conditionalPanel("input.regionradio == 'TRUE'",
                                 sliderInput(
                                   "distenh",
                                   label = strong("Select distance
                                              threshold for merging enhancers"),
                                   min = 0,
                                   max = 3000,
                                   value = 1500
                                   ),
                                 sliderInput(
                                   "distprom",
                                   label = strong("Select distance
                                              threshold for merging promoters"),
                                   min = 0,
                                   max = 3000,
                                   value = 1000
                                   )
                                 ),
                hr(),
                actionButton("buttonannot", strong("Combine and Annotate")),
                hr()
                ),
              HTML("</div>"),

              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              infoBoxOutput("statusbox3", width = NULL),
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
                   400px !important;'>"),
              box(
                width = NULL,
                title = strong("Retrieve Read Counts"),
                h5("This step does the following: "),
                tags$ul(
                  tags$li("Counts the number of reads in each regulatory region
                   for each sample type."),
                  tags$li(" Outputs a denisty plot of the lengths of genomic regions.")
                ),
                hr(),
                h4(" Select Parameters"),
                radioButtons(
                  "chromradio",
                  label = strong("Restrict Analysis to a Single Chromosome?"),
                  choices = list("FALSE" = "FALSE","TRUE" = "TRUE"),
                  selected = "FALSE"
                ),
                conditionalPanel("input.chromradio == 'TRUE'",
                                 uiOutput("chooseChrom")
                ),
                hr(),
                uiOutput("chooseref"),
                hr(),
                actionButton("buttoncounts", strong("Retrieve Counts")),
                hr()
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              infoBoxOutput("statusbox4", width = NULL),

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
                   400px !important;'>"),
              box(
                width = NULL,
                title = strong("Define Altered Regions"),
                h5("This step does the following: "),
                tags$ul(
                  tags$li("Determines which regulatory regions are differentialy
                   altered between sample types.")
                ),
                hr(),
                h4(" Select function parameters used by DESeq2 to calculate
                   adjusted p-values"),

                sliderInput(
                  "alpha",
                  label = strong("Select pvalue cutoff"),
                  min = 0,
                  max = 1,
                  value = 0.01
                ),
                sliderInput(
                  "lfcThreshold",
                  label = strong("Select log2fold change cutoff"),
                  min = 0,
                  max = 5,
                  value = 1,
                  step = 0.1
                ),
                hr(),
                actionButton("buttondefine", strong("Define Altered Regions")),
                hr()
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   500px !important;'>"),
              infoBoxOutput("statusbox5", width = NULL),
              HTML("</div>")
              )
            ),
    tabItem(tabName = "cataltered",
            fluidRow(
              HTML("<div class='col-sm-4' style='min-width:
                   400px !important;'>"),
              box(
                width = NULL,
                title = strong("Categorize Regions"),
                h5("This step does the following: "),
                tags$ul(
                  tags$li("Categorizes altered regulatory regions as experiment-specific,
                reference-specific, or shared."),
                  tags$li(" Outputs a volcano plot highlighting potential
                          altered regulatory elements.")
                ),
                hr(),
                h4("Select parameters that define cell-type specific regulatory regions"),
                sliderInput(
                  "lfcSpecific",
                  label = strong("Select log2fold change cutoff for specific
                             enhancers/promoters"),
                  min = 0,
                  max = 5,
                  value = 1.5,
                  step = 0.1
                ),
                sliderInput(
                  "pvalueSpecific",
                  label = strong("Select pvalue cutoff for specific
                             enhancers/promoters"),
                  min = 0,
                  max = 1,
                  value = 0.01
                ),
                hr(),
                h4("Select parameters that define shared regulatory regions"),
                sliderInput(
                  "lfcShared",
                  label = strong("Select log2fold change cutoff for shared
                             enhancers/promoters"),
                  min = 0,
                  max = 5,
                  value = 1.2,
                  step = 0.1
                ),
                sliderInput(
                  "pvalueShared",
                  label = strong("Select pvalue cutoff for shared enhancers/promoters"),
                  min = 0,
                  max = 1,
                  value = 0.05
                ),
                hr(),
                actionButton("buttoncat", strong("Categorize Altered Regions")),
                hr()
                #,
                # hr(),
                # conditionalPanel("input.buttoncat > 0",
                #                  downloadButton("downloadData",
                #                               strong("Download Track")
                #                               )
                #                  )
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              infoBoxOutput("statusbox5b", width = NULL),
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
                   400px !important;'>"),
              box(
                width = NULL,
                title = strong("Compare Methods"),
                h5("This step does the following: "),
                tags$ul(
                  tags$li("Compares two methods of identifying altered regulatory
                          regions. The first method uses peak presence and
                          associated intensity (i.e chromatin accessibility).
                          The second method uses peak presence only as determined
                          by the peak caller."),
                  tags$li(" Outputs a venn plot")
                ),
                hr(),
                actionButton("buttoncompare", strong("Compare Methods")),
                hr(),
                dataTableOutput("table5")
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              infoBoxOutput("statusbox8", width = NULL),
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
        title = strong("Pathway Enrichment Analysis"),
        width = 12,
        tabPanel("Pathway Enrichment for Molecular Function",
                 fluidRow(
                   infoBoxOutput("statusbox6", width = 11)
                   ),
                 fluidRow(
                   box(
                     title = strong("Pathway Enrichment for Molecular Function"),
                     width = 11,
                     h5("This step does the following: "),
                     tags$ul(
                       tags$li("Determines which pathways are overrepresented in
                             altered promoters and enhancers as returned by
                             GO Enrichment Analysis restricted to the Molecular
                             Function sub-ontology."),
                       tags$li("Outputs a heatmap plot of the enrichment analysis'
                             results.")
                     ),
                     hr(),
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
                   title = strong("Pathway Enrichment for Biological Process"),
                   width = 11,
                   h5("This step does the following: "),
                   tags$ul(
                     tags$li("Determines which pathways are overrepresented in
                             altered promoters and enhancers as returned by
                             GO Enrichment Analysis restricted to the Biological
                             Process sub-ontology."),
                     tags$li("Outputs a heatmap plot of the enrichment analysis'
                             results.")
                   ),
                   hr(),
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
