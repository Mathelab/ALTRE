# App ui

headerbar <- dashboardHeader(
  title = "ALTRE Workflow",
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
      icon = icon("object-ungroup"),
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
      tabName = "greatpathways",
      badgeLabel = "step 8",
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
                dataTableOutput("table2"),
                hr()
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
              hr(),
              HTML("<div class='col-sm-3' style='min-width:
                   350px !important;'>"),
              box(
                title = "Customize Plot",
                width = NULL,
                solidHeader = TRUE,
                uiOutput("choosePalette1"),
                textInput(
                  "consPlotTitle",
                  "Change the main title of the plot",
                  "Peak Counts by Cell Type"
                ),
                hr(),
                textInput(
                  "consPlotylabel",
                  "Change the y-axis label of the plot",
                  "Peak Counts"
                ),
                hr()
              ),
              HTML("</div>"),
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
                h4("Be sure to scroll all the way down for plotting options"),
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
                ),  # end box
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
              hr(),
#              HTML("</div>"),
              HTML("<div class='col-sm-3' style='min-width:
                   350px !important;'>"),
              box(
                title = "Customize Plot",
                width = NULL,
                solidHeader = TRUE,
                uiOutput("choosePalette2"),
                hr(),
                textInput(
                  "combLeftPlotTitle",
                  "Change the main title of the left-hand plot",
                  "Number of REs"
                ),
                hr(),
                textInput(
                  "combRightPlotTitle",
                  "Change the main title of the right-hand plot",
                  "Mean Length of REs"
                ),
                hr(),
                textInput(
                  "combLeftylabel",
                  "Change the y-label of the left-hand plot",
                  "Number of REs"
                ),
                hr(),
                textInput(
                  "combRightylabel",
                  "Change the y-label of the right-hand plot",
                  "Mean Length of REs"
                ),
                hr()
              ),
              HTML("</div>"),
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
                uiOutput("chooseref")
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              box(
                width = NULL,
                #title = "Density Plot",
                highcharter::highchartOutput('densityplot')
              ),
              HTML("<div class='col-sm-3' style='min-width:
                   350px !important;'>"),
              box(
                title = "Customize Plot",
                width = NULL,
                solidHeader = TRUE,
                uiOutput("choosePalette3"),
                hr(),
                textInput(
                  "countsPlotTitle",
                  "Change the main title of the plot",
                  "Density of log2 read counts (normalized by library and region sizes)"
                ),
                hr(),
                textInput(
                  "countsPlotxlabel",
                  "Change the x-axis label of the plot",
                  "Log2 Normalized Read Counts"
                ),
                hr(),
                textInput(
                  "countsPlotylabel",
                  "Change the y-axis label of the plot",
                  "Density"
                ),
                hr()
              ),
              HTML("</div>"),
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
                h5("Note that the total number of categorized peaks may be less than the
                  the total number of peaks evaluated.  This discrepancy is due to DESeq2
                  independent filtering to remove regions with low power.  See DESeq2 documentation
                  for more details."),
                hr(),
                actionButton("buttoncat", strong("Categorize Altered Regions")),
                HTML('<span style="color: #FD3335">(Please Wait! Plot takes ~1 minute to render!)</span>'),
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
                )
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
              HTML("<div class='col-sm-3' style='min-width:
                   350px !important;'>"),
              box(
                title = "Customize Volcano Plot",
                width = NULL,
                solidHeader = TRUE,
                uiOutput("choosePalette4")
              ),
              HTML("</div>"),
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
                  tags$li(" Outputs a set of pie charts.")
                ),
                hr(),
                actionButton("buttoncompare", strong("Compare Methods")),
                hr(),
                dataTableOutput("table5"),
                conditionalPanel("input.buttoncompare > 0",
                                 downloadButton("downloadCompareDT",
                                                strong("Download Data Table")
                                 )
                )
              ),
              box(
                title = "Customize Pie Plot",
                width = NULL,
                solidHeader = TRUE,
                uiOutput("choosePalette5"),
                hr(),
                textInput(
                  "title11",
                  "Change the title of the first plot",
                  "TSS-proximal/Intensity"
                ),
                hr(),
                textInput(
                  "title12",
                  "Change the title of the second plot",
                  "TSS-distal/Intensity"
                ),
                hr(),
                textInput(
                  "title13",
                  "Change the title of the third plot",
                  "All/Intensity"
                ),
                hr(),
                textInput(
                  "title21",
                  "Change the title of the fourth plot",
                  "TSS-proximal/Peak"
                ),
                hr(),
                textInput(
                  "title22",
                  "Change the title of the fifth plot",
                  "TSS-distal/Peak"
                ),
                hr(),
                textInput(
                  "title23",
                  "Change the title of the sixth plot",
                  "All/Peak"
                ),
                hr()
              ),
              HTML("</div>"),
              HTML("<div class='col-sm-7' style='min-width:
                   550px !important;'>"),
              box(
                width = NULL,
                htmlOutput('pieplot')
              ),
              infoBoxOutput("statusbox7", width = NULL),
              HTML("</div>")
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
			      tags$li("Perform pathway analysis with GREAT"),
		              tags$li("Determines which pathways are overrepresented in
                                altered TSS-proximal/distal regions"),
                              tags$li("Outputs a heatmap of the top enriched pathways")
		            ),
			         "Note: You must be connected to the internet for this step.
			           Please be patient. It can take 3-5 minutes to run. Click",
			          tags$a(href = "http://bejerano.stanford.edu/great/public/html/",
			             "here"),
			          "for more information on GREAT.",
		            hr(),
		            actionButton("buttongreat", strong("Run GREAT")),
    			      conditionalPanel("input.buttongreat > 0",
    			                       hr(),
    			                       downloadButton("downloadGREAT",
    			                                      strong("Download Data")
    			                       ))
		          ),
		          HTML("</div>"),
		          HTML("<div class='col-sm-8' style='min-width:
                   750px !important;'>"),
		          box(
		            width = NULL,
		            dataTableOutput("table6")
		          ),
			      box(
               width = NULL,
               #title = "Heat Plot",
               highcharter::highchartOutput('heatplotGREAT')
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
