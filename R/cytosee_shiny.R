#' Opens \code{SC3} results in an interactive session in a web browser.
#'
#' Runs interactive \code{shiny} session of \code{SC3} based on precomputed clusterings.
#'
#' @param object an object of \code{SCESet} class
#'
#' @return Opens a browser window with an interactive \code{shiny} app and visualize
#' all precomputed clusterings.
#'
#' @name cytosee_go
#' @aliases cytosee_go
#' @import shiny
#' @import ggplot2
#' @import visNetwork
#' @import vegan
#' @import pheatmap
#' @import recharts
#' @importFrom  DT renderDataTable dataTableOutput datatable
#' @importFrom shinyjs useShinyjs showElement hideElement alert
#' @importFrom flowCore read.FCS logTransform arcsinhTransform exprs
#' @importFrom shinydashboard dashboardPage sidebarMenu  dashboardHeader dashboardSidebar menuItem box tabBox tabItems tabItem updateTabItems
#' @export

cytosee_go <- function(object){
################ ui ###################
 ui = shinydashboard::dashboardPage(
  skin ="green",
  shinydashboard::dashboardHeader(title = "CytoSEE"),
  shinydashboard::dashboardSidebar(

  #use shinyjs
  useShinyjs(),
  # options
  tags$script(
    HTML(
      "Shiny.addCustomMessageHandler('unbind-DT', function(id) {
        Shiny.unbindAll($('#'+id).find('table').DataTable().table().node());
        })"
    )
  ),
  ##### UI side step one #####
  div(
   id="Step1_sidebar",class="Step1",
   h4("Step.1 Loading and Preview ",style="font-family:Times New Roman"),
   br(),
   div(class="InputFile",
    br(),
    hr(),
    textInput(inputId = "Project_name",value = "CytoSEE",label = "Name the project here:"),
    br(),
    fileInput(inputId = "fcs",
              label = "Only .fcs is available",
              accept=c('.fcs')),
    br()
   ),
  div(
    class="Preview",
    style="display:none",
    uiOutput("XY_val"),
    br(),
    uiOutput("Plot_bin"),
    br(),
    uiOutput("data_transform")
   )
 ),

 ##### UI side Step two #####
 div(id = "Step2_sidebar",class="Step2",
     style="display:none",
     h4("Step.1 Cell clustering ",style="font-family:Times New Roman"),
     radioButtons(inputId = "pre_reduce_method",label = "choose dimension reduce method:",
                 choices=c("ALL","t-SNE","PCA","DiffusionMap","LargeVis"),selected="PCA")
  ),

 ##### UI side Step three #####
 div(id = "Step3_sidebar",class="Step3",
     style="display:none",
     h4("Step 3 Result and Visualization"),
     hr(),
     shinydashboard::sidebarMenu(
       id = "result_display",
       shinydashboard::menuItem("ClusterLabel", icon=icon("tags",lib="font-awesome"), tabName = "ClusterLabel",selected = TRUE),
       shinydashboard::menuItem("ReportTable", icon=icon("sticky-note",lib="font-awesome"), tabName = "ReportTable"),
       shinydashboard::menuItem("Visualization", icon=icon("television",lib="font-awesome"),tabName ="Visualization")
     )
    )
 ),

##### UI Body #####
shinydashboard::dashboardBody(
 ##### UI Body step one #####
 div(
  div(
   class="Introduction",
   h3("Introduction"),
   hr(),
   fluidRow(
    column(
      width=7,
      p("
          Over the past half century, flow cytometry has established itself as the workhorse for high-throughput quantitative analysis of cells and other particles,
          making it one of the earliest omics technologies in retrospect.And a host of new data- analysis tools have emerged,
          creating workflows for processing, analyzing and visualizing complex FCM datasets,
          many of the tools do the similar job but based on quite different algorithms,
          and we cannot determine which is the best because it quite depends on a situation.
          Owning to Comprehensive comparative studies by the Flow Cytometry: Critical Assessment of Population Identification Methods (FlowCAP),
          this project has shown that many of these tools have reached a level of maturity that matches,
          or even surpasses.And FlowCAP indicates using several algorithms in combination might achieve better results. However,
          these workflow or tools have gone mostly unnoticed by non-bioinformaticians especially the immunologists.
          The most probable cause is large amount tools or workflows do not offer friendly UI and require coding experience.
          And many people give it up at the beginning because of the complex installation.
          Another reason may  be the high-throughput nature of flow cytometry,
          combined with the increasing capacity to measure more cell parameters at once,
          is making it more compute-intensive and time-consuming.
          A compute node or at least a workstation is needed to deal with such large dataset.
             "),
      p("
        Here, we present CytoPRO, an integrated R package for pre-processing, clustering,
        visualizing and auto labeling FCS dataset, both work in flow cytometry and Mass cytometry.
        Alongside with a user friendly interface
             ")
        ),
      column(
       width=5,
       plotOutput("")
      )
    )
   ),
   div(
    class="Introduction",
    h3("Some Feature"),
    hr(),
    fluidRow(
    column(
      width=7,
      h5("1. UI based and data driven manipulation:"),
      p("
             CytoPRO begins with a shiny UI application, once the interface launched,
             users can do all steps by mouse clicking without any coding input.
             The FCS data transformed into a cytoSet(an S4 class dataset which contains the raw data and keeps the manipulate records,
             CytoPRO\'s build-in function), via such dataset, A flowchart which contains all the finished analysis can be plotted.
             "),
           h5("2. Multiple and interactive visualization methods:"),
           p("We use the t-SNE`, which is a wildly used algorithm for reducing dimensions and visualize method,
             and also integrate the state-of-art method called largeVis`,  MST, SOM, diffusion map. 2/3D scatter plot. "),
           h5("3. Consensus clustering:"),
           p("As shown in FlowCAP project, using several algorithms in combination might achieve better results."),
           h5("4. Automatics cell population labeling:"),
           p("A new module was developed for automated labeling the clustering result mainly base on flowType and flowCL.
             We reduce the input markers by evaluate the contribution of each markers in clustering,
             then we assess each marker to choose the most reliable phenotype from a series of phenotype in each clusters.
             Finally we get the highest scored labels from the flowCL"),
           h5("5. Other Features"),
           p("elbow method to estimate the K, f-measure, ARI(Adjusted Rand Index), silhouette method to measure the situation of clusters.
             Dynamic report. Parallel computing")
       ),
      column(
        width=5,
        eChartOutput("Statisic_pie",height=400)
     )
    )
   )
  ),
 div(id="Step1_body",
     style="display:none",
     class="Step1",
     fluidRow(
      class="Preview",
      style="display:none",
      column(
        width=12,
        id="lgs",
        box(
         width = 12,
         title = "",
         fluidRow(
          column(
                  width=6,
                  h3("PreViewData"),
                  br(),
                  plotOutput("hex_plot",height = 600),
                  br(),
                  uiOutput("preprocess")
          ),
          column(
                 width=6,
                 h3("Channel Choosing"),
                 br(),
                 box(width = 12,
                 collapsible = FALSE,
                 solidHeader = TRUE,
                 title = "",
                 column(
                    width=12,
                    box(
                      width = 12,
                      DT::dataTableOutput("file")
                    ),
                    actionButton("all_choice","All",width = "25%",style="margin-bottom:5%;margin-left:17%"),
                    actionButton("Done","Done",width = "25%",style="margin-bottom:5%;margin-left:15%;background:#31353e;color:#eeeeee")
                 )
                )
               ),
          # flowAI UI #
          column(
                width=12,
                class="PreProcess",
                style="display:none;",
                tabBox(
                 width = 12,
                 side = "right",
                 title = "",
                 tabPanel(
                  title="flowAI",
                  width=12,
                  fluidPage(
                   box(
                    id="AI_setting",
                    width = 3,
                    title = "AI_setting",
                    h4("Summary:"),
                    textOutput("summaryText1"),
                    textOutput("summaryText2"),
                    hr(),
                    h4("Parameters:"),
                    numericInput("timeLenth", label = h5("Time step (sec)"), value = 0.1, step = 0.1),
                    uiOutput("signalBinSize"),
                    hr(),
                    h4("Download Output:"),
                    br(),
                    downloadButton('downloadQCFCS',   'FCS with QC',style="margin-top:8px;color:#eeeeee;background:#27313e;width:120px"),
                    br(),
                    downloadButton('downloadGoodFCS', 'High Q FCS',style="margin-top:8px;color:#eeeeee;background:#27313e;width:120px"),
                    br(),
                    downloadButton('downloadBadFCS',  'Low Q FCS',style="margin-top:8px;color:#eeeeee;background:#27313e;width:120px")
                  ),
                  tabBox(
                    id="AI_result",width = 9,
                    tabPanel(
                     "Flow Rate",
                     fluidPage(
                      hr(),
                      textOutput("flowRateSummary"),
                      hr(),
                      plotOutput("flowRatePlot",height = 605),
                      hr(),
                      fluidRow(
                       column(4, offset = 1,uiOutput("timeSlider")),
                       column(4, offset = 2,uiOutput("rateSlider"))
                      )
                     )
                    ),
                    tabPanel(
                      "Signal Acquisition",
                      fluidPage(
                       hr(),
                       textOutput("flowSignalSummary"),
                       hr(),
                       uiOutput("signalBinSlider"),
                       hr(),
                       plotOutput("flowSignalPlot", height = "800px")
                       )
                      ),
                    tabPanel(
                      "Dynamic Range",
                      fluidPage(
                       hr(),
                       fluidRow(
                        column(5,textOutput("flowMarginSummary")),
                        column(3, offset = 2,checkboxInput("checkbox", label = "Apply Margins Check", value = TRUE))
                       ),
                       hr(),
                       plotOutput("flowMarginPlot")
                      )
                     )
                    )
                   )
                  )
                 )
                )
              )
             )
            )
           )
  ),
 ##### UI body Step two  #####

 div(id="Step2_body",style="display:none",class="Step2",
          tabsetPanel(
           tabPanel(
            title="DO FAST",
            icon=icon("desktop",lib="font-awesome"),
            fluidRow(
             column(
               width = 4,
               br(),
               box(
                width=12,
                height=600,
                title="flowSOM",
                h4("Introduction"),
                p("FlowSOM, which analyzes Flow or mass cytometry data using a Self-Organizing
                   Map. Using a two-level clustering and star charts, algorithm helps to obtain a clear
                   overview of how all markers are behaving on all cells, and to detect subsets that might be
                   missed otherwise. "),
                sliderInput(inputId = "flowSOM_K",label = "K value means how many groups will be divided.",min = 2,max = 150,step = 1,value = 10,animate = TRUE),
                actionButton(inputId = "flowSOM_run",label = "Calculate",style='background:#31353e;color:#eeeeee;margin-left:25%;width:50%')
                )
               ),
             column(
              width = 4,
              br(),
              box(
                width=12,
                height=600
              )
            ),
             column(
              width = 4,
              br(),
              box(
               width=12,
               height=600
              )
             )
            )
           ),
          tabPanel(
           title="Large-Scale",
           icon=icon("server",lib="font-awesome"),
           fluidRow(
             column(
               width = 4,
               br(),
               box(
                 width=12,
                 height=600,
                 title="DensityCut",
                 h4("Introduction"),
                 p("FlowSOM, which analyzes Flow or mass cytometry data using a Self-Organizing
                   Map. Using a two-level clustering and star charts, algorithm helps to obtain a clear
                   overview of how all markers are behaving on all cells, and to detect subsets that might be
                   missed otherwise. "),
                 sliderInput(inputId = "DensityCut_K",label = "The number of neighbours in building the Knn graph",min = 3,max = 50,step = 1,value = 10,animate = TRUE),
                 actionButton(inputId = "DensityCut_run",label = "Calculate",style='background:#31353e;color:#eeeeee;margin-left:25%;width:50%')
               )
             ),
             column(
               width = 4,
               br(),
               box(
                 width=12,
                 height=600,
                 title="Consboost",
                 h4("Introduction"),
                 p("Consboost, which analyzes Flow or mass cytometry data using a Self-Organizing
                   Map. Using a two-level clustering and star charts, algorithm helps to obtain a clear
                   overview of how all markers are behaving on all cells, and to detect subsets that might be
                   missed otherwise. "),
                 sliderInput(inputId = "Consboost_K",label = "K value means how many groups will be divided.",min = 2,max = 150,step = 1,value = 10,animate = TRUE),
                 actionButton(inputId = "Consboost_run",label = "Calculate",style='background:#31353e;color:#eeeeee;margin-left:25%;width:50%')
               )
             ),
             column(
               width = 4,
               br(),
               box(
                 width=12,
                 height=600,
                 title="Phenograph",
                 h4("Introduction"),
                 p("Phenograph, which analyzes Flow or mass cytometry data using a Self-Organizing
                   Map. Using a two-level clustering and star charts, algorithm helps to obtain a clear
                   overview of how all markers are behaving on all cells, and to detect subsets that might be
                   missed otherwise. "),
                 sliderInput(inputId = "phenograph_K",label = "K value means number of nearest neighbours.",min = 3,max = 50,step = 1,value = 30,animate = TRUE),
                 actionButton(inputId = "phenograph_run",label = "Calculate",style='background:#31353e;color:#eeeeee;margin-left:25%;width:50%')
               )
             )
             )

          ),
          tabPanel(
           title="Small-Scale",
           icon=icon("laptop",lib="font-awesome"),
           fluidRow(
             column(
               width = 4,
               br(),
               box(
                 width=12,
                 height=600,
                 title="flowMeans",
                 h4("Introduction"),
                 p("FlowMeans, which analyzes Flow or mass cytometry data using a Self-Organizing
                   Map. Using a two-level clustering and star charts, algorithm helps to obtain a clear
                   overview of how all markers are behaving on all cells, and to detect subsets that might be
                   missed otherwise. "),
                 sliderInput(inputId = "flowMeans_K",label = "K value means how many groups will be divided.",min = 2,max = 150,step = 1,value = 10,animate = TRUE),
                 sliderInput(inputId = "flowMeans_maxK",label = "maxK value means the Maximun clusters calculated.",min = 2,max = 150,step = 1,value = 20,animate = TRUE),
                 actionButton(inputId = "flowMeans_run",label = "Calculate",style='background:#31353e;color:#eeeeee;margin-left:25%;width:50%')
               )
             )
           )

          )
         )

 ),

 ##### UI body Step three #####
 div(
  id="Step3_body",
  style="display:none",
  class="Step3",
  tabItems(
    tabItem(
      tabName="ClusterLabel",
      div(
        id="change_marker_name",
        h3("Cluster Label Generation",style="color:white"),
        style="background:#dddddd;padding:10px;border-radius:5px;",
        fluidRow(
          column(
           width=4,
           box(
            width=12,
            h4("Note:"),
            p("It is an Cell type identification method for labeling the clusters divided by our clustering methods."),
            p("This approach labeling the clusters based on the relative abundence of cell surface markers e.g.CD38, CCR7."),
            p("Users can change the marker name in the following textbox:"),
            sliderInput(inputId = "CL_MaxHitsPht",label = "MaxHitsPht",min = 1,max = 15,step = 1,value = 5),
            uiOutput("alias")
           )
         ),
         column(
           width=8,
           box(
             width=12,
             height=720,
             id="MarkerLine",
             eChartOutput("Marker_line",height=700)
           ),
           actionButton("clear","clear",width="25%",style="margin-bottom:5%;margin-left:17%"),
           actionButton("label_button",label = "Label",width="25%",style="margin-bottom:5%;margin-left:15%;background:#31353e;color:#eeeeee")
         )
        )
      ),
      div(
        id="celltype_tree_div",
        h3("Clusters Table",style="color:white"),
        style="display:none;background:#dddddd;padding:10px;border-radius:5px;",
        fluidRow(
          column(
            width = 5,
            box(
              width=12,
              DT::dataTableOutput(outputId = "Clusters_table")
              )
          ),
          column(
            width = 7,
            box(
              width = 12,
              height = 720,
              visNetwork::visNetworkOutput("CellType_tree",height = 700,width = "100%")
            ),
            actionButton("Skip",label="Skip",width="25%",style="margin-bottom:5%;margin-left:17%"),
            actionButton("Next",label = "Next",width="25%",style="margin-bottom:5%;margin-left:15%;background:#31353e;color:#eeeeee")
          )
        )
      )
    ),
    tabItem(
     tabName ="ReportTable",
     div(style="background:#dddddd;padding:10px;border-radius:5px;",
      h4("Project Report:"),
      DT::dataTableOutput(outputId = "Project_report",width = "98%"),
      br(),
      h4("Clusters Report:"),
      fluidRow(
        column(
          width = 6,
          DT::dataTableOutput(outputId = "Cluster_report",width = "95%")
        ),
        column(
          width = 6,
          uiOutput(outputId = "confirm_marker_choose"),
          plotOutput(outputId = "cell_type_confirm_plot",height = 330,width = "95%")
        )
      ),
      br(),
      actionButton("previous",label="Previous",width="15%",style="margin-bottom:1%;margin-left:30%"),
      actionButton("confirm",label="Confirm",width="15%",style="margin-bottom:1%;margin-left:10%;background:#31353e;color:#eeeeee")
     )
    ),

    ##### Visualization #####
    tabItem(
     tabName ="Visualization",
     tabBox(
       width=12,
       side = "right",
       id = "DP",
       title = "Visual Display",
       height = "680px",
       selected = "Scatter Plot",

       ### scatter Plot ###
       tabPanel(
         title="Scatter Plot",
         fluidRow(
           column(
             width=9,
             plotOutput(outputId = "scatter_plot",height = 800,width = "100%")
           ),
           column(
             width=3,
             box(
               width = 12,
               title = div(icon("cog"),"\t","Setting"),
               radioButtons(
                 inputId="red_choice",
                 label = "Dimensionality reduction methods:",
                 choices = c("t_SNE","LargeVis","PCA","DiffusionMap"),selected = "PCA",
                 inline = TRUE),
               radioButtons(
                 inputId="Cluster_Marker",
                 label = "Choose a way to display:",
                 choices = c("Clusters","Markers"),selected = "Clusters",
                 inline = TRUE),
               uiOutput("Marker_or_Cluster")
             )
           )
         )
       ),

      ### MST ###
       tabPanel(
         id="MST",
         title = "MST",
         fluidRow(
           column(
             width=8,
             visNetworkOutput("MST",height = 670)
           ),
           column(
             width=4,
             box(
               width = 12,
               uiOutput("MST_panel",height=100)
             )
           )
         )
       ),

       ### Marker heatmap ###
       tabPanel(
         id="Marker_heatmap",
         title="HeatMap",
         fluidRow(
          column(
            width = 8,
            plotOutput("Marker_heatmap",width = "100%",height = 750)
          ),
          column(
            width = 4
          )
         )
       ),

       ### DMT ###
       tabPanel(
         id="DMT",
         title="DMT",
         fluidRow(
          column(
            width = 8,
            visNetwork::visNetworkOutput("DMT",width = "100%",height = 750)
          ),
          column(
            width = 4
          )
         )
       ),
      tabPanel(
        id="Markerline",
        title="Markerline",
        fluidRow(
          column(
            width = 12,
            eChartOutput("cell_type_marker_line",width = "100%",height = 750)
          )
        )
      )
     )
    )
   )
  )
 ) #dashbody end
)  #UI end

######################## server ############################

server = function(input,output,session){
 options(shiny.maxRequestSize=250*1024^2)
 options(rgl.useNULL = TRUE)
 ##### step one #####

 ##### Pie Plot #####
 output$Statisic_pie<-renderEChart({
   Cell_parents = CL_lib[which(CL_lib[,4]=="native cell"),2]
   cell_cont=function(parent){
     if(length(which(CL_lib[,4]==parent))>0){
       children=CL_lib[which(CL_lib[,4]==parent),2]
       count=length(children)
       for(i in 1:length(children)){
         count=count+cell_cont(as.character(children[i]))
       }
       return(count)
     }
     else{
       return(0)
     }
   }
   Count=c()
   cell_type=c()
   for(i in 1:length(Cell_parents)){
     Count=c(cell_cont(as.character(Cell_parents[i]))+1,Count)
     cell_type=c(as.character(Cell_parents[i]),cell_type)
   }
   Cell=data.frame(cell_type,Count)
   Cell=Cell[which(Cell$Count>10),]
   chart = echartr(Cell,cell_type,y=Count,type = "pie") %>%
     setTitle("Cell Ontology main cell types") %>%
     setLegend(show = TRUE,pos = 9) %>%
     setToolbox(language = "en") %>%
     setTheme('roma')
   chart
 })



 # vector for saving choosen markers #
 dy_markers=reactiveValues(data=NULL)

 # a vector for saving file #
 file=reactiveValues(data=NULL)

 # file loading setting #
 observeEvent(input$fcs,{
   showElement("Step1_body",animType ="slide",anim = TRUE,time = 0.5)
   showElement("Step1_bar",animType ="slide",anim = TRUE,time = 0.5)
   hideElement(selector = ".Introduction")
   showElement(selector = ".Preview")
   hideElement(selector = ".InputFile")
   # saving result as cytoSet
   inFile <- input$fcs
   file$data<-read.FCS(inFile$datapath)
   cyto<<-initiflow(fcs.data = exprs(file$data),projectname = input$Project_name)
   cyto@autoLabel<<-FALSE
   cyto@preprocess<<-"None"
 })

 # plot bin set
 output$Plot_bin<-renderUI({
   if(is.null(file$data)){
     return()
   }
   tagList(
     sliderInput(inputId = "bin_val",
                 label = "bin value",
                 min=0,
                 max=512,
                 value = 256,
                 step=16),
     br(),
     sliderInput(inputId = "scale_val",
                 label = "scale value",
                 min=10,
                 max=500,
                 value = 30,
                 step=10)
   )
 })

 #data Preview
 output$hex_plot<-renderPlot({
   if(is.null(file$data)){
     return()
   }
   if(is.null(input$transform)){
     return()
   }
   data<- file$data
   if (input$transform=="log10"){
     tran <- logTransform(transformationId="log10-transformation", logbase=10, r=1, d=1)
     tranlist_log <- transformList(colnames(data), tran)
     data<-transform(data, tranlist_log)
     cyto@transform_method<<-"log10"
   }

   else if (input$transform=="arcsinh"){
     tran <- arcsinhTransform("arcsinh-transformation",a = 1,b=1,c=0)
     tranlist_arcsinh <- transformList(colnames(data), tran)
     data<-transform(data, tranlist_arcsinh)
     cyto@transform_method<<-"arcsinh"
   }
   else{
     cyto@transform_method<<-"None"
   }
   data_trans$data<-as.data.frame(exprs(data))
   scale=input$scale_val
   cyto@fcs.data<<-data_trans$data
   data=data_trans$data
   ggplot(data,
          aes(x=data[input$plotx],
              y=data[input$ploty]))+
     geom_hex(bins=input$bin_val)+
     theme_bw()+
     theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
     xlab(input$plotx)+
     ylab(input$ploty)+
     theme(panel.border  = element_blank())+
     theme(axis.line = element_line(colour = "black"))+
     scale_fill_gradientn(colours = c("black","darkred","orange","yellow","white"),
                          limits=c(0,scale),
                          name="Count of cells\n")
 })

 # file description #
 output$file <- DT::renderDataTable({
   inFile <- file$data
   if (is.null(inFile))
     return(NULL)
   fcs=file$data
   # define the base view table for files
   table<-data.frame(as.data.frame(fcs@parameters@data$name),
                     as.data.frame(fcs@parameters@data$desc),
                     as.data.frame(fcs@parameters@data$range),
                     as.data.frame(fcs@parameters@data$minRange),
                     as.data.frame(fcs@parameters@data$maxRange),
                     row.names = NULL)

   colnames(table)<-c("name","desc","range","minRange","maxRange")

   DT::datatable(
     table,
     class = "hover",
     style = "bootstrap",
     options = list(dom="rpit",paging=FALSE,scrollY=450),
     selection = list(mode="multiple",
                      selected=r_all_choice$data,
                      target="row")
   )
 })

 # choose the display markers
 output$XY_val <- renderUI({
   fcs<-file$data
   if (is.null(fcs)){
     return()
   }
   tagList(
     h3("Plot Setting"),
     selectInput("plotx", "plot-x",
                 colnames(fcs),
                 selected = colnames(fcs)[1]
     ),
     selectInput("ploty", "plot-y",
                 colnames(fcs),
                 selected = colnames(fcs)[2]
     )
   )
 })

 # data_transform
 data_trans<-reactiveValues(data=NULL)
 output$data_transform <- renderUI({
   if(is.null(file$data)){
     return()
   }
   tagList(
     radioButtons(inputId = "transform",label = "data transform",
                  choices = c("log10","arcsinh","None"),
                  selected = "None",inline = TRUE)
   )
 })

 # function of "all" button #
 r_all_choice=reactiveValues(data=NULL)
 observeEvent(input$all_choice,{
   r_all_choice$data=NULL
   r_all_choice$data=c(1:length(file$data@parameters@data$name))
 })


 # function of "Done" button #
 observeEvent(input$Done,{
   inFile<-input$fcs
   if (is.null(inFile)){
     return("A .fcs is needed")
   }
   fcs=read.FCS(inFile$datapath)
   channel_id = input$file_rows_selected
   dy_markers$data=colnames(exprs(fcs))[as.integer(channel_id)]
   if(length(dy_markers$data)==0){
     alert("No Channel was chosen, at least one Channel is required!")
     return()
   }
   hideElement(selector = ".InputFile")
   showElement(selector = ".Preview")
   hideElement(selector = ".Step1",animType ="fade",anim = TRUE,time = 0.3)
   hideElement("Step1_bar")
   showElement(selector = ".Step2",animType ="slide",anim = TRUE,time = 1)
   showElement("Step2_bar")
   cyto@channel.use<<-channel_id
 })

 # preprocess UI #
 output$preprocess<-renderUI({
   data=data_trans$data
   if(is.null(data)){
     return()
   }
   tagList(
     h3("PreProcess"),
     checkboxInput(inputId = "Choose_flowAI",label = "flowAI",value = FALSE),
     p("(flowAI is a quality control tools which can be efficent in figuring out outliers and eliminate them.)")
   )
 })

 ####### flowAI Servers #######

 ## use flowAI selection
  observeEvent(input$Choose_flowAI,{
    if(input$Choose_flowAI=="TRUE"){
      showElement(selector = ".PreProcess",anim = TRUE,animType = "slide")
      #hideElement("lgs",anim = TRUE,animType = "slide")
      cyto@fcs.data<<-exprs(checkRes()[[3]])
      cyto@preprocess<<-"flowAI"
    }
    else{
      hideElement(selector = ".PreProcess",anim = TRUE,animType = "slide")
      cyto@fcs.data<<-data_trans$data
    }
  })

  ## load flowset data
  set <- reactive({
    isolate({fcsFiles <- input$fcs
    if (is.null(fcsFiles))
      return(NULL)
    set <- read.FCS(fcsFiles$datapath)
    set@description$FILENAME <- fcsFiles$name})
    return(set)
  })

  ## time channel name
  timeChannel <- reactive({
    if(is.null(set()))
      return(NULL)
    x <- set()
    time <- findTimeChannel(x)
    return(time)
  })

  ## time step
  timeStep <- reactive({
    if(is.null(set()))
      return(NULL)
    word <- which(grepl("TIMESTEP", names(set()@description),
                        ignore.case = TRUE))
    timestep <- as.numeric(set()@description[[word[1]]])
    if( !length(timestep) ){
      warning("The timestep keyword was not found in the FCS file and it was set to 0.01. Graphs labels indicating time might not be correct", call. =FALSE)
      timestep <- 0.01
    }
    return(timestep)
  })


  TimeChCheck <- reactive({
    if (!is.null(timeChannel())) {
      if (length(unique(exprs(set())[, timeChannel()])) == 1){
        TimeChCheck <- "single_value"
      }else{
        TimeChCheck <- NULL
      }
    }else{
      TimeChCheck <- "NoTime"
    }
    return(TimeChCheck)
  })


  ## order fcs expression according acquisition time
  ordFCS <- reactive({
    if(is.null(set()))
      return(NULL)
    if(is.null(TimeChCheck())){
      ordFCS <- ord_fcs_time(set(), timeChannel())
    }else{
      ordFCS <- set()
    }
    return(ordFCS)
  })


  ## signal bin size UI
  output$signalBinSize <- renderUI({
    if(is.null(set())){
      optSize <- NULL
      maxSize <- Inf
    }else{
      maxSize <- nrow(ordFCS())
      optSize <- min(max(1, floor(maxSize/100)), 500)
    }
    numericInput("signalBinSize", label = h5("Number of events per bin:"),
                 value = optSize, min = 1, max = maxSize)
  })


  ## cell quality check
  cellCheck <- reactive({
    if(is.null(ordFCS()))
      return(NULL)
    if(is.null(TimeChCheck())){
      flowRateData <- flow_rate_bin(ordFCS(), second_fraction = input$timeLenth,
                                    timeCh = timeChannel(), timestep = timeStep())
    }else{
      flowRateData <- list()
    }
    flowSignalData <- flow_signal_bin(ordFCS(), channels = NULL,
                                      binSize = input$signalBinSize, timeCh = timeChannel(),
                                      timestep = timeStep(), TimeChCheck = TimeChCheck() )

    flowMarginData <- flow_margin_check(ordFCS())

    res <- list(flowRateData, flowSignalData, flowMarginData)
    return(res)
  })


  ## flow rate time slider UI and check sliders. if they are null, a default value is returned for the QC
  sliders <- reactive({
    flowRateData <- cellCheck()[[1]]
    flowSignalData <- cellCheck()[[2]]
    return(c(
      min(flowRateData$frequencies[,3]) - 0.1,
      max(flowRateData$frequencies[,3]) + 0.1,
      min(flowRateData$frequencies[,4]) - 10,
      max(flowRateData$frequencies[,4]) + 10,
      0,
      nrow(flowSignalData$exprsBin) + 1)
    )
  })

  output$timeSlider <- renderUI({
    if(is.null(set()) || is.null(cellCheck()) || !is.null(TimeChCheck()))
      return(NULL)
    sliderInput("timeSlider", strong("Time cut:"),
                min = sliders()[1], max = sliders()[2],
                value = c(sliders()[1], sliders()[2]), step = 0.1)
  })
  timeSlider <- reactive({
    if(is.null(input$timeSlider)){
      return(c(sliders()[1], sliders()[2]))
    }else{
      return(c(input$timeSlider[1],  input$timeSlider[2]))
    }

  })

  output$rateSlider <- renderUI({
    if(is.null(set()) || is.null(cellCheck()) || !is.null(TimeChCheck()))
      return(NULL)
    sliderInput("rateSlider", strong("Flow rate cut:"),
                min = sliders()[3], max = sliders()[4],
                value = c(sliders()[3], sliders()[4]), step = 0.1)
  })
  rateSlider <- reactive({
    if(is.null(input$rateSlider)){
      flowRateData <- cellCheck()[[1]]
      return(c(sliders()[3], sliders()[4]))
    }else{
      return(c(input$rateSlider[1],  input$rateSlider[2]))
    }

  })

  output$signalBinSlider <- renderUI({
    if(is.null(set()) || is.null(cellCheck()))
      return(NULL)
    sliderInput("signalBinSlider", strong("Signal acquisition cut:"), width = "90%",
                min = sliders()[5], max = sliders()[6],
                value = c(sliders()[5], sliders()[6]), step = 1)
  })
  signalSlider <- reactive({
    if(is.null(input$signalBinSlider)){
      return(c(sliders()[5], sliders()[6]))
    }else{
      return(c(input$signalBinSlider[1],  input$signalBinSlider[2]))
    }
  })


  ## plot
  output$flowRatePlot <- renderPlot({
    if(is.null(ordFCS()) || is.null(cellCheck()) || !is.null(TimeChCheck()))
      return(NULL)
    flowRateData <- cellCheck()[[1]]
    frp <- flow_rate_plot(flowRateData, input$rateSlider[1], input$rateSlider[2],
                          input$timeSlider[1], input$timeSlider[2])
    print(frp)
  })

  output$flowSignalPlot <- renderPlot({
    if(is.null(set()) || is.null(cellCheck()))
      return(NULL)
    flowSignalData <- cellCheck()[[2]]
    fsp <- flow_signal_plot(flowSignalData, input$signalBinSlider[1], input$signalBinSlider[2])
    print(fsp)
  })

  output$flowMarginPlot <- renderPlot({
    if(is.null(set()) || is.null(cellCheck()))
      return(NULL)
    flowMarginData <- cellCheck()[[3]]
    fmp <- flow_margin_plot(flowMarginData, input$signalBinSize)
    print(fmp)
  })



  ## check results
  checkRes <- reactive({
    if(is.null(set()) || is.null(cellCheck()))
      return(NULL)

    ordFCS <- ordFCS()
    totalCellNum <- nrow(ordFCS)
    origin_cellIDs <- 1:totalCellNum
    if(is.null(TimeChCheck())){
      FlowRateQC <- flow_rate_check(cellCheck()[[1]], rateSlider()[1], rateSlider()[2],
                                    timeSlider()[1], timeSlider()[2])
    }else{
      FlowRateQC <- origin_cellIDs
    }
    FlowSignalQC <- flow_signal_check(cellCheck()[[2]], signalSlider()[1], signalSlider()[2])

    if(input$checkbox[1] == TRUE){
      FlowMarginQC <- cellCheck()[[3]]$goodCellIDs
    }else{
      FlowMarginQC <- origin_cellIDs
    }

    goodCellIDs <- intersect(FlowRateQC, intersect(FlowSignalQC, FlowMarginQC))
    badCellIDs <- setdiff(origin_cellIDs, goodCellIDs)

    flowRatePerc <- 1 - length(FlowRateQC)/length(origin_cellIDs)
    flowSignalPerc <- 1 - length(FlowSignalQC)/length(origin_cellIDs)
    flowMarginPerc <- 1 - length(FlowMarginQC)/length(origin_cellIDs)
    totalBadPerc <- length(badCellIDs)/length(origin_cellIDs)

    params <- parameters(ordFCS)
    keyval <- keyword(ordFCS)
    sub_exprs <- exprs(ordFCS)

    good_sub_exprs <- sub_exprs[goodCellIDs, ]
    goodfcs <- flowFrame(exprs = good_sub_exprs, parameters = params, description = keyval)

    bad_sub_exprs <- sub_exprs[badCellIDs, ]
    badfcs <- flowFrame(exprs = bad_sub_exprs, parameters = params,description = keyval)

    tempQCvector <- cellCheck()[[2]]
    QCvector <- tempQCvector$cellBinID[,"binID"]
    QCvector[badCellIDs] <- runif(length(badCellIDs), min=10000, max=20000)
    QCfcs <- addQC(QCvector, sub_exprs, params, keyval)

    return(list(totalCellNum, totalBadPerc, goodfcs, badfcs,
                flowRatePerc, flowSignalPerc, flowMarginPerc, QCfcs))
  })

  ## summary text
  output$summaryText1 <- renderText({
    if(is.null(checkRes()))
      return(NULL)
    paste0("Total number of events: ", checkRes()[[1]])
  })

  output$summaryText2 <- renderText({
    if(is.null(checkRes()))
      return(NULL)
    paste0("Percentage of low-Q events: ", round(checkRes()[[2]]*100,2), "%")
  })

  output$flowRateSummary <- renderText({
    if(is.null(checkRes()))
      return(NULL)
    if(is.null(TimeChCheck())){
      paste0("Percentage of low-Q events in flow rate check: ", round(checkRes()[[5]]*100,2), "%")
    }else if(!is.null(TimeChCheck()) && TimeChCheck() == "NoTime"){
      "It is not possible to recreate the flow rate because the time channel is missing."
    }else if(!is.null(TimeChCheck()) && TimeChCheck() == "single_value"){
      "It is not possible to recreate the flow rate because the time channel contains a single value."
    }
  })

  output$flowSignalSummary <- renderText({
    if(is.null(checkRes()))
      return(NULL)
    paste0("Percentage of low-Q events in signal acquisition check: ", round(checkRes()[[6]]*100,2), "%")
  })

  output$flowMarginSummary <- renderText({
    if(is.null(checkRes()))
      return(NULL)
    paste0("Percentage of low-Q events in dynamic range check: ", round(checkRes()[[7]]*100,2), "%")
  })

  file_base <- reactive({
    file_ext <- description(ordFCS())$FILENAME
    file_base <- sub("^([^.]*).*", "\\1", file_ext)
    return(file_base)
  })

  ## download processed FCS files
  output$downloadGoodFCS <- downloadHandler(
    filename = function(){
      paste0(file_base(), "_HighQ.fcs")
    },

    content = function(file){
      data <- checkRes()[[3]]
      if(is.null(data)){
        return(NULL)
      }
      write.FCS(data, file)
      #tar(tarfile = file, files = tempdir)
    }
  )

  output$downloadBadFCS <- downloadHandler(
    filename = function(){
      paste0(file_base(), "_LowQ.fcs")
    },

    content = function(file){
      data <- checkRes()[[4]]
      if(is.null(data)){
        return(NULL)
      }
      write.FCS(data, file)
      #tar(tarfile = file, files = tempdir)
    }
  )

  ## download processed FCS files
  output$downloadQCFCS <- downloadHandler(
    filename = function(){
      paste(file_base(), "_QC.fcs", sep='')
    },

    content = function(file){
      data <- checkRes()[[8]]
      if(is.null(data)){
        return(NULL)
      }
      write.FCS(data, file)
      #tar(tarfile = file, files = tempdir)
    })

  ##### Step 2 Cell clustering #####

  ### flowSOM ###
  observeEvent(input$flowSOM_run,{
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                  file=flowFrame(exprs = as.matrix(cyto@fcs.data))
                  incProgress(1/10,message="Run flowSOM...")
                  ClusterID=cytosee_flowSOM(File = file,K = input$flowSOM_K,colsToUse = cyto@channel.use)
                  ClusterID=as.data.frame(ClusterID)
                  data=cbind(ClusterID,cyto@fcs.data[,cyto@channel.use])
                  cyto@label<<-ClusterID
                  Cluster_id <- unique(cyto@label)
                  cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                  cyto@ClusterID<<-ClusterID
                  cyto@clust_method<<-"flowSOM"
                  incProgress(1/10,message="Run DMT...")
                  cyto@dmt<<-DMT(data)

                  incProgress(1/10,message="Deduce the Dimension...")
                  cyto@dim.red<<-pro_reduce_dim(cyto@fcs.data[,cyto@channel.use],method=input$pre_reduce_method)
                  incProgress(7/10,message="Complete!")
                 })
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })

  ### DensityCut ###
  observeEvent(input$DensityCut_run,{
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   file= cyto@fcs.data[,cyto@channel.use]
                   incProgress(1/10,message="Run DesnityCut")
                   ClusterID=cytosee_densitycut(file,NumC = input$DensityCut_K)
                   ClusterID=as.data.frame(ClusterID)
                   data=cbind(ClusterID,cyto@fcs.data[,cyto@channel.use])
                   cyto@label<<-ClusterID
                   Cluster_id <- unique(cyto@label)
                   cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                   cyto@ClusterID<<-ClusterID
                   cyto@clust_method<<-"DensityCut"
                   incProgress(1/10,message="Run DMT...")
                   cyto@dmt<<-DMT(data)

                   incProgress(1/10,message="Reduce the Dimension...")
                   cyto@dim.red<<-pro_reduce_dim(cyto@fcs.data[,cyto@channel.use],method=input$pre_reduce_method)
                   incProgress(7/10,message="Complete!")
                 })
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })


  ### Consboost ###
  observeEvent(input$Consboost_run,{
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   file= cyto@fcs.data[,cyto@channel.use]
                   incProgress(1/10,message="Run Consboost")
                   ClusterID=cytosee_consboost(file,maxK=input$Consboost_K,K=input$Consboost_K)
                   ClusterID=as.data.frame(ClusterID)
                   data=cbind(ClusterID,cyto@fcs.data[,cyto@channel.use])
                   cyto@label<<-ClusterID
                   Cluster_id <- unique(cyto@label)
                   cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                   cyto@ClusterID<<-ClusterID
                   cyto@clust_method<<-"Consboost"
                   incProgress(1/10,message="Run DMT...")
                   cyto@dmt<<-DMT(data)

                   incProgress(1/10,message="Reduce the Dimension...")
                   cyto@dim.red<<-pro_reduce_dim(cyto@fcs.data[,cyto@channel.use],method=input$pre_reduce_method)
                   incProgress(7/10,message="Complete!")
                 })
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })


  ### FlowMeans ###
  observeEvent(input$flowMeans_run,{
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   file= cyto@fcs.data[,cyto@channel.use]
                   incProgress(1/10,message="Run FlowMeans")
                   ClusterID=cytosee_flowMeans(file,MaxN=input$flowMeans_maxK,NumC=input$flowMeans_K)
                   ClusterID=as.data.frame(ClusterID)
                   data=cbind(ClusterID,cyto@fcs.data[,cyto@channel.use])
                   cyto@label<<-ClusterID
                   Cluster_id <- unique(cyto@label)
                   cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                   cyto@ClusterID<<-ClusterID
                   cyto@clust_method<<-"FlowMeans"
                   incProgress(1/10,message="Run DMT...")
                   cyto@dmt<<-DMT(data)

                   incProgress(1/10,message="Reduce the Dimension...")
                   cyto@dim.red<<-pro_reduce_dim(cyto@fcs.data[,cyto@channel.use],method=input$pre_reduce_method)
                   incProgress(7/10,message="Complete!")
                 })
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })


  ### Phenograph ###
  observeEvent(input$phenograph_run,{
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   file= cyto@fcs.data[,cyto@channel.use]
                   incProgress(1/10,message="Run Phenograph")
                   ClusterID=cytosee_phenograph(file,K=input$phenograph_K)
                   ClusterID=as.data.frame(ClusterID)
                   data=cbind(ClusterID,cyto@fcs.data[,cyto@channel.use])
                   cyto@label<<-ClusterID
                   Cluster_id <- unique(cyto@label)
                   cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                   cyto@ClusterID<<-ClusterID
                   cyto@clust_method<<-"Phenograph"
                   incProgress(1/10,message="Run DMT...")
                   cyto@dmt<<-DMT(data)

                   incProgress(1/10,message="Reduce the Dimension...")
                   cyto@dim.red<<-pro_reduce_dim(cyto@fcs.data[,cyto@channel.use],method=input$pre_reduce_method)
                   incProgress(7/10,message="Complete!")
                 })
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })

  ##### Step 3 Labeling and Visualization #####

  ### Label the clusters###
  output$alias<-renderUI({
    if(is.null(cyto@dmt)){
      return()
    }
    markers=c()
    for(i in c(1:length(cyto@dmt)-1)){
      re=unlist(strsplit(cyto@dmt[[i+1]]$label,"[<||>]"))
      markers[i]=re[1] #in order to delet the root
    }
    markers=unique(markers)
    taglist=list()
    j=1

    for(i in c(1:length(markers))){
      marker_id=sprintf("M_%s",j)
      tag=fluidRow(
        column(
          width=6,
          textInput(inputId =marker_id,
                    label = markers[i],
                    value = r_marker_value$data[i],
                    width = 150)
        )
      )
      taglist=c(taglist,list(tag))
      j=j+1
      #eval(parse(text=sprintf("select_marker=c(input$%s,select_marker)",i)))
    }
    div(
      style='overflow-y:auto;
      overflow-x:hidden;
      height:550px',
      taglist
    )
  })

  ##### label button #####
  observeEvent(input$label_button,{
    if(is.null(cyto@dmt)){
      return()
    }
    markers=c()
    for(i in c(1:length(cyto@dmt)-1)){
      re=unlist(strsplit(cyto@dmt[[i+1]]$label,"[<||>]"))
      markers[i]=re[1] #in order to delet the root
    }
    markers=unique(markers)
    select_marker=list()
    j=1
    for(i in markers){
      marker_id=sprintf("M_%s",j)
      if(eval(parse(text = sprintf("!is.na(input$%s)",marker_id)))){
        eval(parse(text=sprintf("select_marker=c(select_marker,input$%s)",marker_id)))
      }
      else{
        eval(parse(text=sprintf("select_marker=c(select_marker,%s)",NA)))
      }
      j=j+1
    }

    ### data transfom for labeling ###
    marker2=unlist(markers)
    if(is.null(select_marker)){
      alert("Input area shouldn't be empty!")
      return()
    }
    else{
      manual_des=unlist(select_marker)
    }

    marker_change=c()
    for(i in c(1:length(marker2))){
      marker_change[marker2[i]]=manual_des[i]
    }

    Cluster_marker=c()
    Cluster_m2c=c()
    cnt=0
    for(i in c(1:length(cyto@dmt))){
      if(length(cyto@dmt[[i]]$clus)==1){
        tmp=cyto@dmt[[i]]
        str=""
        while(tmp$id!=0){
          str=paste0(tmp$label,",",str)
          tmp=cyto@dmt[[tmp$previd+1]]
        }
        cnt=cnt+1
        Cluster_marker[cnt]=str
        Cluster_m2c[cnt]=cyto@dmt[[i]]$clus
      }
    }
    for(i in c(1:length(marker_change))){
      Cluster_marker=gsub(pattern = names(marker_change[i]),replacement = marker_change[i],Cluster_marker,fixed = TRUE)
    }
    re=gsub(pattern = " ",replacement = "",Cluster_marker)
    re=gsub(pattern = "<.*?,",replacement = "<",re)
    re=gsub(pattern = ">.*?,",replacement = ">",re)

    for(s in c(1:length(marker_change))){
      word=marker_change[s]
      for(i in c(1:length(re))){
        store=re[i]
        plus=0
        minus=0

        while(length(grep(pattern = paste0(word,"<"),x = store,perl = FALSE))!=0 || length(grep(pattern = paste0(word,">"),x = store,perl = FALSE))!=0 ){
          word_p=paste0(word,">")
          word_m=paste0(word,"<")

          if(length(grep(word_p,store))!=0){
            store=sub(pattern=word_p,replacement="",store,perl = FALSE)
            plus=plus+1
          }
          if(length(grep(word_m,store))!=0){
            store=sub(pattern = word_m,replacement="",store,perl = FALSE)
            minus=minus+1
          }
        }
        if(plus==0&&minus==0){
          next()
        }
        if(plus-minus>=2){
          re[i]=paste0(store,word,"++")
        }
        else if(plus-minus==1){
          re[i]=paste0(store,word,"+")
        }
        else if(plus==minus){
          re[i]=paste0(store,word,"+")
        }
        else if(plus-minus==-1){
          re[i]=paste0(store,word,"-")
        }
        else{
          re[i]=paste0(store,word,"--")
        }
      }
    }
    CL_label=list()
    withProgress(
      message = "Run LocCL",
      value=0,
      {
    for(i in 1:length(re)){
      incProgress(1/length(re))
      result=LocCL(MarkerList = re[i],MaxHitsPht = input$CL_MaxHitsPht)
      ClusterID=Cluster_m2c[i]
      CL_label[[i]]=list("result"=result,"ClusterID"=ClusterID)
    }
   }
  )
    cyto@CL_label<<-CL_label
    cyto@autoLabel<<-TRUE
    hideElement(id="change_marker_name")
    showElement(id="celltype_tree_div")
  })

  ##### clear button #####
  r_marker_value=reactiveValues(data=c())
  observeEvent(input$clear,{
    if(is.null(cyto@dmt)){
      return()
    }
    markers=c()
    for(i in c(1:length(cyto@dmt)-1)){
      markers[i]=unlist(strsplit(cyto@dmt[[i+1]]$label,"[<||>]"))[1] #in order to delet the root
    }
    markers=unique(markers)
    for(i in c(1:length(markers))){
      r_marker_value$data[i]="1"
      r_marker_value$data[i]=""
    }
  })

  ##### Marker_line #####
  output$Marker_line<-renderEChart({
    cl_average_vec<-getAvg(fcsData=cyto@fcs.data[,cyto@channel.use],ClusterID=cyto@label)
    ClusterID<-c()
    Channels<-c()
    Value<-c()
    for(i in 1:length(cl_average_vec)){
      for(j in 1:length(cl_average_vec[,1])){
        ClusterID<-c(row.names(cl_average_vec[j,]),ClusterID)
        Channels<-c(colnames(cl_average_vec)[i],Channels)
        Value<-c(cl_average_vec[j,i],Value)
      }
    }
    Line<-data.frame(ClusterID,Channels,Value)
    echartr(Line,x=Channels,y = Value,series = ClusterID,type = "curve") %>%
      setLegend(selected = ClusterID[1],pos=6) %>%
      setTheme('macarons2') %>%
      setToolbox(language = "en") %>%
      setTitle(title = "Marker expression level in each clusters",pos = 12)
  })


  # helper function for making checkbox in CL label datatable
  shinyInput_CL = function(FUN, len, id, ...){
    inputs = character(len)
    for (i in seq_len(len)) {
      inputs[i] = as.character(FUN(paste0(id, i),paste0(cyto@CL_label[[i]]$result$Labels),...))
    }
    inputs
  }

  # helper function for making checkbox in Cluster report datatable
  shinyInput_Cluster = function(FUN, len, id,Label, ...){
    inputs = character(len)
    for (i in seq_len(len)) {
      inputs[i] = as.character(FUN(paste0(id, i),paste0(Label[i]),...))
    }
    inputs
  }


  # helper function for reading checkbox in datatable
  shinyValue = function(id, len) {
    unlist(lapply(seq_len(len), function(i) {
      value = input[[paste0(id, i)]]
      if (is.null(value))
        NA
      else
        value
    }))
  }

  ##### Clusters_tabel #####
  output$Clusters_table<-DT::renderDataTable({
    Cluster_id = unique(cyto@label)
    size = c()
    for(i in 1:length(Cluster_id[,1])){
      size = c(size,length(cyto@label[which(cyto@label==Cluster_id[i,]),]))
    }
    table=data.frame(Cluster_id,size,"Label Choose"=shinyInput_CL(
      selectInput,
      nrow(Cluster_id),
      "selecter_",
      width = "390px",
      label = NULL
    )
   )
    data.frame(table)
  }, selection = list(mode='single',selected=1), server = FALSE, escape = FALSE,rownames = FALSE, options = list(
    dom = "rpit",
    paging = FALSE,
    scrollY=650,
    preDrawCallback=JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
    drawCallback=JS("function() { Shiny.bindAll(this.api().table().node()); } ")
  ))

  ##### Cell type tree builder #####
  output$CellType_tree <-renderVisNetwork({
    id=input$Clusters_table_rows_selected
    if(is.null(id)){
      return("you can choose a cluster to display")
    }
    Nodes<-cyto@CL_label[[id]]$result$Nodes
    Links<-cyto@CL_label[[id]]$result$Links
    visNetwork(Nodes, Links) %>%
      visEdges(arrows = "from",color="#cccccc") %>%
      visNodes(size=15) %>%
      visOptions(highlightNearest = list(enabled=TRUE,algorithm="hierarchical",
                                         degree=list(from=0,to=50))) %>%
      visHierarchicalLayout(treeSpacing=300,nodeSpacing=150,
                            levelSeparation = 100,direction = "DU",
                            blockShifting = FALSE,sortMethod = "directed")
  })

  ##### Skip #####
  observeEvent(input$Skip,{
    Cluster_id <- unique(cyto@label)
    cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
    updateTabItems(session,"result_display",selected="ReportTable")
  })


  ##### Next button #####
  observeEvent(input$Next,{
    Cluster_id <- unique(cyto@label)
    cyto@ID2CL <<- data.frame(ID=Cluster_id,CD_Label = shinyValue("selecter_", nrow(Cluster_id)))
    updateTabItems(session,"result_display",selected="ReportTable")
  })


  ##### Project_report #####
  output$Project_report<-DT::renderDataTable({
    projectName=cyto@projectname
    cell_num<-length(cyto@fcs.data[,1])
    cluster_num<-length(unique(cyto@label)[,1])
    transform_method<-cyto@transform_method
    channels<-paste(colnames(cyto@fcs.data)[cyto@channel.use],collapse = ",\t")
    clust_method<-cyto@clust_method
    preprocess<-cyto@preprocess
    AutoLabel<-cyto@autoLabel
    table<-data.frame(c("ProjectName","Cell num","num of clusters","Channels","Transform method","Clust Method","Preprocess Method","AutoLabel"),
                      c(projectName,cell_num,cluster_num,channels,transform_method,clust_method,preprocess,AutoLabel)
    )
    colnames(table)=c("Characters","Values")

    DT::datatable(
      table,
      class = "hover",
      style = "default",
      rownames=FALSE,
      escape = FALSE,
      selection="none",
      options = list(dom="rpt",paging = FALSE,ordering=FALSE,searching=FALSE)
    )
  })

  ##### Cluster_report #####
  output$Cluster_report<-DT::renderDataTable({
    input$Next
    Cluster_id = unique(cyto@ClusterID)
    ID2CL=cyto@ID2CL
    inputLabel = c()
    for(i in 1:length(Cluster_id[,1])){
      inputLabel =c(inputLabel,as.character(ID2CL[which(ID2CL[,1]==Cluster_id[i,1]),2]))
    }
    size = c()
    for(i in 1:length(Cluster_id[,1])){
      size = c(size,length(cyto@label[which(cyto@label==Cluster_id[i,]),]))
    }
    table=data.frame(Cluster_id,size,"Label Input"=shinyInput_Cluster(
      textInput,
      len=nrow(Cluster_id),
      id="text_",
      Label = inputLabel,
      width = "395px",
      label = NULL
    )
   )
    data.frame(table)
  }, selection = list(selected=1,mode="single"), server = FALSE, escape = FALSE,rownames = FALSE, options = list(
    dom = "rpit",
    paging = FALSE,
    scrollY=330,
    preDrawCallback=JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
    drawCallback=JS("function() { Shiny.bindAll(this.api().table().node()); } ")
  )
 )

  ##### cell_type_confirm_plot #####
  output$cell_type_confirm_plot<-renderPlot({
    clusters=unlist(cyto@label)
    unique_cluster=unique(clusters)
    clusters[which(clusters!=unique_cluster[input$Cluster_report_rows_selected])]="Others"
    data=cyto@fcs.data[,cyto@channel.use]
    Dimension1<<-data[,input$confirm_marker_x]
    Dimension2<<-data[,input$confirm_marker_y]
    ggplot(data = data,aes(x=Dimension1,y=Dimension2))+geom_point(aes(colour=clusters))+theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      theme(panel.border  = element_blank())+
      theme(axis.line = element_line(colour = "black"))+
      xlab(input$confirm_marker_x)+
      ylab(input$confirm_marker_y)+
      scale_colour_manual(values = c("#3dccfb","#c2c2c2"))
  })

  ##### cell_type_confirm_plot_marker choose #####
  output$confirm_marker_choose<-renderUI({
    fluidRow(
      column(
        width = 6,
        selectInput(inputId = "confirm_marker_x",
                    label = "X",
                    selected = colnames(cyto@fcs.data)[cyto@channel.use[1]],
                    choices = colnames(cyto@fcs.data)[cyto@channel.use])
             ),
      column(
        width = 6,
        selectInput(inputId = "confirm_marker_y",
                    label = "Y",
                    selected = colnames(cyto@fcs.data)[cyto@channel.use[2]],
                    choices = colnames(cyto@fcs.data)[cyto@channel.use])
      )
    )
  })

  ##### confirm button #####
  observeEvent(input$confirm,{
    Cluster_id = unique(cyto@ClusterID)
    cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = shinyValue("text_", nrow(Cluster_id)))
    ID2CL<-cyto@ID2CL
    label=cyto@label
    for(i in 1:length(Cluster_id[,1])){
      label[which(label[,1]==Cluster_id[i,1]),1]<-as.character(ID2CL[which(ID2CL[,1]==Cluster_id[i,1]),2])
    }
    cyto@label<<-label
    updateTabItems(session,"result_display",selected="Visualization")
  })

  ##### previous button #####
  observeEvent(input$previous,{
    updateTabItems(session,"result_display",selected="ClusterLabel")
  })

  ##### DMT plot#####
  output$DMT<-renderVisNetwork({
    if(is.null(cyto@dmt)){
      return()
    }

    #vector for nodes
    id=c(0)
    label=c("Root")
    title=c("Root")
    size=c(length(cyto@fcs.data[,1]))
    group=c("Root")

    #vector for links
    from=c()
    to=c()
    for(i in c(2:length(cyto@dmt))){
      if(length(cyto@dmt[i][[1]]$clus)>1){
        group=c("Others",group)
        title=c(paste0("ClusterID:","Others"),title)
      }
      else{
        group=c("Clusters",group)
        title=c(paste0("ClusterID:",cyto@dmt[i][[1]]$clus[1]),title)
      }
      label=c(cyto@dmt[i][[1]]$label,label)
      id=c(cyto@dmt[i][[1]]$id,id)
      size=c(cyto@dmt[i][[1]]$size,size)
      from=c(cyto@dmt[i][[1]]$previd,from)
      to=c(cyto@dmt[i][[1]]$id,to)
    }
    links=data.frame(from=from,to=to,length=2)
    nodes=data.frame(id=id,label=label,size=log10(size)*10+3,title=title,group=group)
    visNetwork(nodes,links, width = "100%",height = "100%") %>%
      visEdges(arrows = "to") %>%
        visOptions(highlightNearest = list(enabled=TRUE,algorithm="hierarchical",
                                           degree=list(from=100,to=0)),
                   selectedBy = list(variable="title",multiple=TRUE)) %>%
        visGroups(groupname = "Root", color ="#6f89ab") %>%
        visGroups(groupname = "Clusters", color = "#f0eccc") %>%
        visGroups(groupname = "Others", color = "#ffdca3")

    })

  ##### cell_type_marker_line #####
  output$cell_type_marker_line<-renderEChart({
    cl_average_vec<-getAvg(fcsData=cyto@fcs.data[,cyto@channel.use],ClusterID=cyto@label)
    ClusterID<-c()
    Channels<-c()
    Value<-c()
    for(i in 1:length(cl_average_vec)){
      for(j in 1:length(cl_average_vec[,1])){
        ClusterID<-c(row.names(cl_average_vec[j,]),ClusterID)
        Channels<-c(colnames(cl_average_vec)[i],Channels)
        Value<-c(cl_average_vec[j,i],Value)
      }
    }
    Line<-data.frame(ClusterID,Channels,Value)
    echartr(Line,x=Channels,y = Value,series = ClusterID,type = "curve") %>%
      setLegend(selected = ClusterID[1],pos=3) %>%
      setTheme('macarons2') %>%
      setToolbox(language = "en") %>%
      setGrid(width = 1100,height = 650) %>%
      setTitle(title = "Marker expression level in each cell types",pos = 12)
  })

 ##### MST plot #####
  output$MST<-renderVisNetwork({
    cl_average_vec<-getAvg(fcsData=cyto@fcs.data[,cyto@channel.use],ClusterID=cyto@label)
    dis <- vegdist(cl_average_vec)
    tr <- spantree(dis)
    Color_group=c("#dcdcdc","#bfe2bc","#9df5b3","#4cecd7","#00deff")
    scale_all=cl_average_vec[input$MST_markers]
    scale_max=max(as.matrix(scale_all))
    scale_min=min(as.matrix(scale_all))
    scale_diff=(scale_max-scale_min)/4
    colortable=c()
    for(i in c(1:length(tr$labels))){
      c=scale_all[which(row.names(scale_all)==tr$labels[i]),]
      for(j in c(1:5)){
        if(c<scale_min+j*scale_diff && c>=scale_min+scale_diff*(j-1)){
          colortable=c(colortable,Color_group[j])
          break()
        }
      }
    }
    node=data.frame(id=1:length(tr$labels),label=tr$labels,color=colortable)
    link=data.frame(from=2:(length(tr$kid)+1),to=tr$kid,color=c("#cccccc"))
    visNetwork(node, link)

  })
 ### MST pannel ###
  output$MST_panel<-renderUI({
    tagList(
      selectInput(inputId = "MST_markers",label = "Markers",choices = as.matrix(colnames(cyto@fcs.data[,cyto@channel.use])),selected = colnames(cyto@fcs.data[,cyto@channel.use])[1])
    )
  })

 ##### Marker heatmap #####
  output$Marker_heatmap<-renderPlot({
    cluster_av=getAvg(cyto@fcs.data[,cyto@channel.use],cyto@label)
    pheatmap::pheatmap(cluster_av,color = colorRampPalette(c("#ffffd9","#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58"))(100))
  })

 ##### scatter plot #####
  output$scatter_plot<-renderPlot({
    if(input$red_choice=="t_SNE"){
      if(!is.null(cyto@dim.red$tsne2d)){
        data=as.data.frame(cyto@dim.red$tsne2d$Y)
      }
      else{
        return(ggplot()+ggtitle("This method is not choosen! Please try others"))
      }
    }
    else if(input$red_choice=="LargeVis"){
      if(!is.null(cyto@dim.red$largeVis)){
        data=as.data.frame(t(cyto@dim.red$largeVis$coords))
      }
      else{
        return(ggplot()+ggtitle("This method is not choosen! Please try others"))
      }
    }
    else if(input$red_choice=="PCA"){
      if(!is.null(cyto@dim.red$PCA)){
        data=as.data.frame(cyto@dim.red$PCA$scores[,c(1,2)])
      }
      else{
        return(ggplot()+ggtitle("This method is not choosen! Please try others"))
      }
    }
    else if(input$red_choice=="DiffusionMap"){
      if(!is.null(cyto@dim.red$diffusion2d)){
        data=as.data.frame(cyto@dim.red$diffusion2d@eigenvectors)
      }
      else{
        return(ggplot()+ggtitle("This method is not choosen! Please try others"))
      }
    }
    Dimension1=data[,1]
    Dimension2=data[,2]
    if(input$Cluster_Marker=="Markers"){
      Expression=cyto@fcs.data[input$Marker_exp_sca_choice]
      ggplot(data = data,aes(x=Dimension1,y=Dimension2))+geom_point(aes(colour=Expression))+theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())+
        theme(panel.border  = element_blank())+
        theme(axis.line = element_line(colour = "black"))+
        scale_color_gradientn(colours = c("#ffffd9","#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58"))
    }
    else if(input$Cluster_Marker=="Clusters"){
      if(input$Cluster_sca_choice=="ALL"){
        clusters=as.factor(unlist(cyto@label))
        ggplot(data = data,aes(x=Dimension1,y=Dimension2))+geom_point(aes(colour=clusters))+theme_bw()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
          theme(panel.border  = element_blank())+
          theme(axis.line = element_line(colour = "black"))
      }
      else{
        clusters=unlist(cyto@label)
        clusters[which(clusters!=input$Cluster_sca_choice)]="Others"
        ggplot(data = data,aes(x=Dimension1,y=Dimension2))+geom_point(aes(colour=clusters))+theme_bw()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
          theme(panel.border  = element_blank())+
          theme(axis.line = element_line(colour = "black"))+
          scale_colour_manual(values = c("#3dccfb","#c2c2c2"))
      }
    }
  })

 ###  markers or clusters event###
 observeEvent(input$Cluster_Marker,{
  if(input$Cluster_Marker=="Markers"){
    showElement(id="Marker_exp_sca")
    hideElement(id="Cluster_sca")
  }
  else if(input$Cluster_Marker=="Clusters"){
    showElement(id="Cluster_sca")
    hideElement(id="Marker_exp_sca")
  }
 })

 ### marker or cluster event UI ###
 output$Marker_or_Cluster=renderUI({
  tagList(
    div(
      id="Cluster_sca",
      selectInput(inputId = "Cluster_sca_choice",
                  label = "you can select a cluster for displaying:",
                  choices = c("ALL",unique(cyto@label)),selected = "ALL")
    ),
    div(
      id="Marker_exp_sca",style="display:none",
      selectInput(inputId = "Marker_exp_sca_choice",
                  label = "you can select a Marker for displaying:",
                  choices = c(colnames(cyto@fcs.data[,cyto@channel.use])),
                  selected=cyto@fcs.data[,cyto@channel.use][1])
    )
   )
 })

 }
   shinyApp(ui=ui, server = server)
}
