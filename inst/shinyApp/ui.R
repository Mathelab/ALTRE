# App ui

headerbar <- dashboardHeader(
  title = h4("ALTRE Analysis Workflow"),
  titleWidth = 270,
  dropdownMenu(
    type = "notifications",
    notificationItem(
      text = "Plots might take some time to display",
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
      "Identify Consensus Peaks",
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
      "Identify Altered Regions",
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
    ),
    menuItem(
      "Find Enriched Pathways",
      icon = icon("gears"),
      tabName = "greatpathways",
      badgeLabel = "step 8a",
      badgeColor = "green"
    ),
    menuItem(
    actionButton("buttonstop", strong("Click to Exit Shiny App")),
    icon = icon("sign-out")
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
                   650px !important;'>"),
              box(
                title = strong("Load Metadata Spreadsheet"),
                width = NULL,
                solidHeader = TRUE,
                tags$b("Please be sure that all files noted in the CSV file,
                         including the CSV file, are in the same folder."),
                h5("This step does the following: "),
                tags$ul(
                  tags$li("Loads a metadata spreadsheet with a CSV file extention."),
		  tags$li("CSV file must contain the following columns:
                        'bamfiles' (name of alignment files),
                        'replicate' (replicate number),
                        'sample' (name of sample),
                        'peakfiles' (name of peak files)"),
                  tags$li("Prints out the contents of the file.")
                ),
                hr(),
	        strong("Load CSV File:"),
		br(),
		shinyFilesButton('file',
		                 'Select File',
		                 'Provide CSV File to Load Data',
		                 FALSE),
  		hr(),
              textOutput("getlocalpath"),
		  dataTableOutput("table1")
                ),
              HTML("</div>"),
              HTML("<div class='col-sm-6' style='min-width:
                   600px !important;'>"),
              infoBoxOutput("statusbox1", width = NULL),
              HTML("</div>")
              )
            ),
    tabItem(tabName = "definerep",
            fluidRow(
              HTML("<div class='col-sm-4' style='min-width:
                   500px !important;'>"),
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
                tags$p(" Note: You can save the plots by clicking on the right mouse button
                  and selecting 'save image as'"),
                hr(),
                actionButton("buttonmerge",
                             strong("Load Files then Merge Replicates")),
                hr(),
                numericInput(
                  "numOverlap",
                  "Choose the Minimum Number of Overlaping Genomic Regions N",
                  2,
                  min = 2,
                  max = 10
                ),
                hr(),
                dataTableOutput("table2")
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              hr(),
              box(
                #title = "Barplot",
                width = NULL,
                solidHeader = TRUE,
                highcharter::highchartOutput("barplot")
              ),
              infoBoxOutput("statusbox2", width = NULL),
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
                          on whether each region is TSS-proximal or TSS-distal
                          as determined by user-defined distance from a
                          transcription start site.")
                ),
                hr(),
                actionButton("buttonannot", strong("Combine and Annotate")),
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
                                                  for merging TSS-proximal and TSS-distal"),
                                   min = 0,
                                   max = 3000,
                                   value = 0
                                   )
                                 ),
                conditionalPanel("input.regionradio == 'TRUE'",
                                 sliderInput(
                                   "distTSSdist",
                                   label = strong("Select distance
                                              threshold for merging TSS-distal regions"),
                                   min = 0,
                                   max = 3000,
                                   value = 1500
                                   ),
                                 sliderInput(
                                   "distTSSprox",
                                   label = strong("Select distance
                                              threshold for merging TSS-proximal regions"),
                                   min = 0,
                                   max = 3000,
                                   value = 1000
                                   )
                                 ),
                conditionalPanel("input.buttonannot> 0",
                                 hr(),
                                 downloadButton("downloadAnnotate",
                                                strong("Download CSV File with
                                                       Annotated Regions ")
                                 )
                ),
                hr()
                ),
              HTML("</div>"),

              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              box(
                #title = "Barplot",
                width = NULL,
                solidHeader = TRUE,
                htmlOutput('annotatebarplot')
              ),
              box(
                #title = "Barplot",
                width = NULL,
                solidHeader = TRUE,
                hr(),
                dataTableOutput("table3")
              ),
              infoBoxOutput("statusbox3", width = NULL),
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
                  tags$li(" Outputs a density plot of the widths of genomic regions.")
                ),
                hr(),
                actionButton("buttoncounts", strong("Retrieve Counts")),
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
                hr()
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              box(
                width = NULL,
                #title = "Density Plot",
                highcharter::highchartOutput('densityplot')
              ),
              infoBoxOutput("statusbox4", width = NULL),
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
                actionButton("buttondefine", strong("Define Altered Regions")),
                hr(),
                h4(" Select function parameters used by DESeq2 to calculate
                   adjusted p-values"),

                sliderInput(
                  "alpha",
                  label = strong("Select p-value cutoff"),
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
                actionButton("buttoncat", strong("Categorize Altered Regions")),
                hr(),
                h4("Select parameters that define cell-type specific regulatory regions"),
                sliderInput(
                  "lfcSpecific",
                  label = strong("Select log2fold change cutoff for specific
                             TSS-proximal/distal regions"),
                  min = 0,
                  max = 5,
                  value = 1.5,
                  step = 0.1
                ),
                sliderInput(
                  "pvalueSpecific",
                  label = strong("Select p-value cutoff for specific
                             TSS-proximal/distal regions"),
                  min = 0,
                  max = 1,
                  value = 0.01
                ),
                hr(),
                h4("Select parameters that define shared regulatory regions"),
                sliderInput(
                  "lfcShared",
                  label = strong("Select log2fold change cutoff for shared
                             TSS-proximal/distal regions"),
                  min = 0,
                  max = 5,
                  value = 1.2,
                  step = 0.1
                ),
                sliderInput(
                  "pvalueShared",
                  label = strong("Select p-value cutoff for shared TSS-proximal/distal regions"),
                  min = 0,
                  max = 1,
                  value = 0.05
                ),
                conditionalPanel("input.buttoncat > 0",
                                 hr(),
                                 downloadButton("downloadBED",
                                                strong("Download BED File")
                                 )
                ),
                hr()
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              box(
                width = NULL,
                #title = "Volcano Plot",
                solidHeader = TRUE,
                htmlOutput('volcano')
              ),
              box(
                width = NULL,
                #title = "Volcano Plot",
                highcharter::highchartOutput('boxplotCounts')
              ),
              box(
                width = NULL,
                #title = "Volcano Plot",
                dataTableOutput("table4")
              ),
              infoBoxOutput("statusbox6", width = NULL),
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
                          associated intensity (i.e. chromatin accessibility).
                          The second method uses peak presence only as determined
                          by the peak caller."),
                  tags$li(" Outputs a Venn diagram.")
                ),
                hr(),
                actionButton("buttoncompare", strong("Compare Methods")),
                hr(),
                dataTableOutput("table5"),
                # conditionalPanel("input.buttoncompare > 0",
                #                  hr(),
                #                  downloadButton("downloadCompareDT",
                #                                 strong("Download Data Table")
                #                  )
                # ),
                hr()
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              box(
                width = NULL,
                htmlOutput('vennplot')
              ),
              infoBoxOutput("statusbox7", width = NULL),
              HTML("</div>")
              )
            ),
    tabItem(
      tabName = "pathways",
      tabBox(
        title = strong("Pathway Enrichment Analysis"),
        width = 12,
        tabPanel("(1) Pathway Enrichment for Molecular Function",
                 fluidRow(
                   HTML("<div class='col-sm-4' style='min-width:
                   400px !important;'>"),
                   box(
                     width = NULL,
                     title = strong("Pathway Enrichment for Molecular Function"),
                     h5("This step does the following: "),
                     tags$ul(
                       tags$li("Determines which pathways are overrepresented in
                             altered TSS-proximal/distal regions as returned by
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
                       label = strong("Select p-value cutoff"),
                       min = 0,
                       max = 1,
                       value = 0.01
                     ),
                     conditionalPanel("input.buttonpathwayMF> 0",
                                      hr(),
                                      downloadButton("downloadPathwayMF",
                                                     strong("Download Zip
                                                            File with MF
                                                            Analysis Results")
                                      )
                     ),
                     hr()
                   ),
                   HTML("</div>"),
                   HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
                   box(
                     width = NULL,
                     #title = "Heat Plot",
                     highcharter::highchartOutput('heatplotMF')
                   ),
                   infoBoxOutput("statusbox8a", width = NULL),
                   HTML("</div>")
                 )
                 ),
       tabPanel("(2) Pathway Enrichment for Biological Process",
                fluidRow(
                  HTML("<div class='col-sm-4' style='min-width:
                       400px !important;'>"),
                  box(
                    width = NULL,
                    title = strong("Pathway Enrichment for Biological Process"),
                    h5("This step does the following: "),
                    tags$ul(
                      tags$li("Determines which pathways are overrepresented in
                             altered TSS-proximal/distal regions as returned by
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
                      label = strong("Select p-value cutoff"),
                      min = 0,
                      max = 1,
                      value = 0.01
                    ),
                    conditionalPanel("input.buttonpathwayBP> 0",
                                     hr(),
                                     downloadButton("downloadPathwayBP",
                                                    strong("Download Zip
                                                            File with BP
                                                            Analysis Results")
                                     )
                    ),
                    hr()
                  ),
                  HTML("</div>"),
                  HTML("<div class='col-sm-7' style='min-width:
                       550px !important;'>"),
                  box(
                    width = NULL,
                    #title = "Heat Plot",
                    highcharter::highchartOutput('heatplotBP')
                  ),
                  infoBoxOutput("statusbox8b", width = NULL),
                  HTML("</div>")
                  )
                ),
       tabPanel("(3) Pathway Enrichment for Cellular Component",
                fluidRow(
                  HTML("<div class='col-sm-4' style='min-width:
                       400px !important;'>"),
                  box(
                    width = NULL,
                    title = strong("Pathway Enrichment for Cellular Component"),
                    h5("This step does the following: "),
                    tags$ul(
                      tags$li("Determines which pathways are overrepresented in
                              altered TSS-proximal/distal regions as returned by
                              GO Enrichment Analysis restricted to the
                              Cellular Component sub-ontology."),
                      tags$li("Outputs a heatmap plot of the enrichment analysis'
                              results.")
                      ),
                    hr(),
                    actionButton("buttonpathwayCC",
                                 strong("Run CC Pathway Enrichment")),
                    hr(),
                    sliderInput(
                      "pathpvaluecutoffCC",
                      label = strong("Select p-value cutoff"),
                      min = 0,
                      max = 1,
                      value = 0.01
                    ),
                    conditionalPanel("input.buttonpathwayCC> 0",
                                     hr(),
                                     downloadButton("downloadPathwayCC",
                                                    strong("Download Zip
                                                            File with CC
                                                            Analysis Results")
                                     )
                    ),
                    hr()
                      ),
                  HTML("</div>"),
                  HTML("<div class='col-sm-7' style='min-width:
                       550px !important;'>"),
                  box(
                    width = NULL,
                    #title = "Heat Plot",
                    highcharter::highchartOutput('heatplotCC')
                  ),
                  infoBoxOutput("statusbox8c", width = NULL),
                  HTML("</div>")
                  )
      )
       )
      ),
		tabItem(tabName = "greatpathways",
		        fluidRow(
		          HTML("<div class='col-sm-3' style='min-width:
		               400px !important;'>"),
		          box(
		            width = NULL,
		            title = strong("GREAT Pathways"),
		            h5("This step does the following: "),
		            tags$ul(
		              tags$li("Runs GREAT."),
		              tags$li(" Outputs a ..")
		            ),
		            hr(),
		            actionButton("buttongreat", strong("Run GREAT")),
		            hr()
		          ),
		          HTML("</div>"),
		          HTML("<div class='col-sm-8' style='min-width:
                   750px !important;'>"),
		          box(
		            width = NULL,
		            dataTableOutput("table6")
		          ),
		          infoBoxOutput("statusbox9", width = NULL),
		          HTML("</div>")
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
