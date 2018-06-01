#' Opens \code{CytoSEE} workflow in an interactive session in a web browser.
#'
#' Runs interactive \code{shiny} session of \code{CytoSEE}
#'
#' @param object an object of \code{cytosee} class
#'
#' @return Opens a browser window with an interactive \code{shiny} app and visualize
#' all precomputed clusterings.
#'
#' @name cytosee_gui
#' @aliases cytosee_gui
#'
#' @import shiny
#' @import visNetwork
#' @import pheatmap
#' @import recharts
#' @import hexbin
#' @import ggplot2
#' @importFrom shinyalert shinyalert useShinyalert
#' @importFrom shinyjs useShinyjs showElement hideElement alert
#' @importFrom flowCore read.FCS logTransform arcsinhTransform exprs transformList transform colnames flowFrame parameters keyword description<- description fsApply
#' @importFrom shinydashboard dashboardPage sidebarMenu  dashboardHeader dashboardSidebar menuItem box tabBox tabItems tabItem updateTabItems
#' @export

cytosee_gui <- function(object){

##===========================shiny UI part=================================

ui = shinydashboard::dashboardPage(
  skin ="blue",
  shinydashboard::dashboardHeader(title = "CytoSEE"),

  ##### sider bar #####
  shinydashboard::dashboardSidebar(

    # use shinyjs and shinyalert
    useShinyjs(),
    useShinyalert(),
    disable = TRUE,

  tags$script(
    HTML(
      "Shiny.addCustomMessageHandler('unbind-DT', function(id) {
        Shiny.unbindAll($('#'+id).find('table').DataTable().table().node());
        })"
    )
  ),
 div(id = "Step3_sidebar",class="Step3",
     style="display:none",
     h4("Step 3 Result and Visualization"),
     hr(),
     shinydashboard::sidebarMenu(
       id = "result_display",
       shinydashboard::menuItem("ClusterLabel", icon=icon("tags",lib="font-awesome"), tabName = "ClusterLabel",selected = TRUE),
       shinydashboard::menuItem("ClusterinfoTable", icon=icon("sticky-note",lib="font-awesome"), tabName = "ClusterinfoTable"),
       shinydashboard::menuItem("Visualization", icon=icon("television",lib="font-awesome"),tabName ="Visualization")
     )
    )
 ),

##### UI Body #####
shinydashboard::dashboardBody(
 ##### UI Body step one #####
 div(
  class="Introduction",
  div(
   hr(),
   fluidRow(
     column(
       width = 12,
       h1("CytoSEE:",br(),span("Shiny based platform for single cells analyzing",style="font-family:'Georgia'; font-size:20px;"),
        br(),
        style =
        "font-family: 'Bell MT';
        font-size:50px;
        color: #000; text-align: center;
        padding: 30px;
        border-radius:5px")
     )
    )
   ),
   div(
    hr(),
    fluidRow(
    column(
      width = 4,
      box(
        title="Single-file analysis",
        status = "primary",
        width = 12,
        solidHeader=TRUE,
        collapsible=TRUE,
        div(class="InputFile",
            br(),
            h4("Note:"),
            h5("This module is an entrance for single file analysis and display."),
            hr(),
            textInput(inputId = "Project_name1",value = "CytoSEE",label = "Name the project here:"),
            br(),
            fileInput(inputId = "fcs",
                      label = "Only .fcs is permitted",
                      accept=c('.fcs')),
            br()
      )
     )
    ),
    column(
      width = 4,
      box(
        title="Rdata display",
        status = "primary",
        width = 12,
        solidHeader=TRUE,
        collapsible=TRUE,
        div(class="InputFile",
            br(),
            h4("Note:"),
            h5("This module is an entrance for dataset which is already analyzed by CytoSEE."),
            hr(),
            textInput(inputId = "Project_name2",value = "CytoSEE",label = "Name the project here:"),
            br(),
            fileInput(inputId = "Rdata",
                      label = "Only .Rdata is permitted",
                      accept=c('.Rdata')),
            br()
        )
      )
    ),
    column(
      width = 4,
      box(
        title="Multi-file analysis",
        status = "primary",
        width = 12,
        solidHeader=TRUE,

        collapsible=TRUE,
        div(class="InputFile",
            br(),
            h4("Note:"),
            h5("This module is an entrance for multi-file analysis."),
            hr(),
            textInput(inputId = "Project_name3",value = "CytoSEE",label = "Name the project here:"),
            br(),
            fileInput(inputId = "zip",
                      label = "Only .zip is permitted",
                      accept=c('.zip')),
            br()
        )
      )
     )
    ),
    h4("Contact:mchen@zju.edu.cn",br(),sprintf("VERSION: %s",packageVersion("cytosee")),
         style =
        "font-family: 'Bell MT';
         font-size:19px;
         color: #000; text-align: center;
         padding: 10px;
         margin-top:5%")
   )
  ),

 div(
     id="Step1_body",
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
         title = "Preview",
         fluidRow(
          column(
                  width=8,
                  br(),
                  fluidRow(
                    column(
                      width = 6,
                      plotOutput("hex_plot",height = 450)
                    ),
                    column(
                      width = 6,
                      plotOutput("X_density",height = 220),
                      plotOutput("Y_density",height = 220)
                    )
                  ),
                  fluidRow(
                    column(
                      width=12,
                      uiOutput("XY_val"),
                      uiOutput("Plot_bin"),
                      actionButton("back_homepage","back",width="20%",style="margin-top:3%;margin-left:20%"),
                      actionButton("previewer_next","Next",width = "20%",style="margin-top:3%;margin-left:20%;background:#31353e;color:#eeeeee")
                    )
                  )
          ),
          column(
                 width=4,
                 br(),
                 div(
                   column(
                     width=12,
                     box(
                       title = "Markers for clustering",
                       width = 12,
                       DT::dataTableOutput("file"),
                       actionButton("all_choice","All",width = "20%",style="margin-top:1%;margin-left:75%")
                     ),
                     box(
                       title = "Markers for data transformation",
                       width = 12,
                       uiOutput("data_transform"),
                       hr(),
                       div(
                        style="overflow:auto;height:160px",
                        helpText("Deselect markers that do not need to be transformed"),
                        uiOutput(outputId = "Choose_transform_channel")
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


 #####============== UI body PreProcess ==============#####
 div(
   id="choose_QC",class="Preprocess",style="display:none",
   column(
     width = 4
   ),
   column(
     width = 4,
     mainPanel(
       width = 12,
       style="background:#fff;border-radius:5px;padding:10px;border-style: solid;border-color:#aaa;border-width:1px",
       radioButtons(inputId = "select_QC",label = "Choose a method for QC",choices = c("flowAI-auto","flowAI-manual"),inline=TRUE),
       helpText("choose a QC method here.if'skip' button is clicked QC step will be skipped"),
       br(),
	     actionButton(inputId = "back_preview",label = "back", width="20%",style="margin-bottom:5%;margin-left:10%"),
       actionButton(inputId = "select_QC_skip",label = "Skip",width = "20%",style="margin-bottom:5%;margin-left:10%"),
       actionButton(inputId = "select_QC_next",label = "Next",width = "20%",style="margin-bottom:5%;margin-left:10%;background:#31353e;color:#eeeeee")
     )
   )
 ),
 # flowAI UI
 div(id="flowAI_UI",class="Preprocess",style="display:none",
   box(
     width = 12,
     column(
       width=12,
       h4("flowAI"),
       br(),
       class="PreProcess",
       fluidPage(
         tabBox(
           id="AI_result",width = 8,
           tabPanel(
             "Flow Rate",
             fluidPage(
               hr(),
               textOutput("flowRateSummary"),
               hr(),
               plotOutput("flowRatePlot"),
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
               div(
                 style="overflow:auto;height:450px",
                 uiOutput("flowSignalPlotUI",width = "98%")
               )
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
         ) ,
        box(
           id="AI_setting",
           width = 4,
           title = "AI_setting",
           h4("Summary:"),
           textOutput("summaryText1"),
           textOutput("summaryText2"),
           hr(),
           h4("Parameters:"),
           numericInput("timeLenth", label = h5("Time step (sec)"), value = 0.1, step = 0.1),
           uiOutput("signalBinSize"),
           br(),
           actionButton(inputId = "AI_Skip",label = "Skip",width = "20%",style="margin-bottom:5%;margin-left:20%"),
           actionButton(inputId = "AI_Confirm",label = "Confirm",width = "20%",style="margin-bottom:5%;margin-left:20%;background:#31353e;color:#eeeeee")
         )
       )
     )
   )
 ),



 ##### UI body Step two  #####

 div(id="Step2_body",style="display:none",class="Step2",
        div(
             column(
               width = 1,
               br()
             ),
             column(
               width = 2,
               uiOutput("File_info")
             ),
             column(
               width = 6,
               mainPanel(
                 width = 12,
                 id = "Step2_sidebar",class="Step2",
                 style="display:none;background:#fff;border-radius:5px;padding:10px;border-style: solid;border-color:#aaa;border-width:1px",
                 h3("Cell clustering",style ="text-align: center;color:#111111"),
                 h3(icon("cog"),"Setting"),
                 column(
                   width = 6,
                   radioButtons(inputId = "pre_reduce_method",
                                label = "Dimensionality reduction method:",
                                choices=c("t-SNE","PCA","LargeVis"),
                                selected="PCA",
                                inline = TRUE),
                   selectInput(inputId = "select_clust_method",
                               label = "Clustering methods:",
                               selected = "consboost",
                               choices = c("Consboost","DensityCut","Flowmeans","FlowSOM","SamSPECTRAL","Phenograph"))
                 ),
                 column(
                   width= 6,
                   radioButtons(inputId = "pre_mst_method",
                                label = "MST methods:",
                                choices = c("flowSOM"),
                                selected = "flowSOM",
                                inline = TRUE),
                   fileInput(inputId = "import_cluster",
                             label = "Custom Clusters: '.csv' file is allowed.",
                             accept = c(".csv"))
                 ),
                 div(
                   div(
                     actionButton(inputId = "back_preprocess",label = "back",width = "15%" ,style="margin-left:80%;margin-bottom:5px")
                   )
                 ),
                 hr(),
                 div(
                   id="FlowSOM_clust",
                   style="display:none",
                   h3(icon("desktop"),"FlowSOM"),
                   column(
                     width = 12,
                     h4(strong("Description:"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif"),
                     p("Method to run general FlowSOM workflow. Will scale the data and uses consensus meta-clustering by default."),
                     h4(strong("speed:"),icon("battery-full"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif")
                   ),
                   br(),
                   column(
                     width = 4,
                     sliderInput(inputId = "flowSOM_K",label = "maxK value",min = 2,max = 150,step = 1,value = 10)
                   ),
                   column(
                     width = 4,
                     textInput(inputId = "flowSOM_nclus",label = "nClus", value = "NULL")
                   ),
                   column(
                     width = 4,
                     selectInput(inputId = "flowSOM_metaclust_method",label = "metaclustering method", choices = c("metaClustering_consensus","metaClustering_hclust","metaClustering_kmeans","metaClustering_som") ,selected = "metaClustering_consensus")
                   ),
                   actionButton(inputId = "flowSOM_run",label = "Calculate",style='margin-left:80%;width:15%;background:#31353e;color:#eeeeee')
                 ),
                 div(
                   id="DensityCut_clust",
                   style="display:none",
                   h3(icon("desktop"),"DensityCut"),
                   column(
                     width = 12,
                     h4(strong("Description:"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif"),
                     p("densityCut first roughly estimates the densities of data points from a K-nearest neighbour graph and then refines the densities via a random walk.
                       A cluster consists of points falling into the basin of attraction of an estimated mode of the underlining density function.
                       A post-processing step merges clusters and generates a hierarchical cluster tree.
                       The number of clusters is selected from the most stable clustering in the hierarchical cluster tree. Experimental results on ten synthetic benchmark datasets and two microarray gene expression datasets demonstrate that densityCut performs better than state-of-the-art algorithms for clustering biological datasets."),
                     h4(strong("speed:"),icon("battery-half"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif")
                   ),
                   br(),
                   column(
                     width=6,
                     numericInput(inputId = "DensityCut_K",label = "Knn K",min = 10,max = 1000,step = 10,value = 30),
                     sliderInput(inputId = "DensityCut_alpha",label = "alpha",min = 0,max = 1,step = 0.1,value = 0.90)
                   ),
                   column(
                     width=6,
                     numericInput(inputId = "DensityCut_maxit",label = "maxit",min = 10,value = 50),
                     sliderInput(inputId = "DensityCut_eps",label = "eps",min = 1e-10,max = 0.05,value = 1e-5)
                   ),
                   actionButton(inputId = "DensityCut_run",label = "Calculate",style='margin-left:80%;width:15%;background:#31353e;color:#eeeeee')
                   ),
                 div(
                   id="Consboost_clust",
                   h3(icon("desktop"),"Consboost"),
                   column(
                     width = 12,
                     h4(strong("Description:"),sytle="font-size:18px;font-family:'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif"),
                     p("adaboost;consensus clustering;down-sampling"),
                     h4(strong("speed:"),icon("battery-quarter"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif")
                   ),
                   br(),
                   column(
                     width = 4,
                     sliderInput(inputId = "Consboost_K",label = "K value",min = 2,max = 150,step = 1,value = 10)
                   ),
                   column(
                     width = 4,
                     sliderInput(inputId = "Consboost_sample_n",label = "sample_n",min = 1000,max = 5000,step = 100,value = 3000)
                   ),
                   column(
                     width = 4,
                     numericInput(inputId = "Consboost_n_core",label = "core",min = 1,value = 3)
                   ),
                   actionButton(inputId = "Consboost_run",label = "Calculate",style='margin-left:80%;width:15%;background:#31353e;color:#eeeeee')
                 ),
                 div(
                   id="Phenograph_clust",
                   style="display:none",
                   h3(icon("desktop"),"Phenograph"),
                   column(
                     width = 12,
                     h4(strong("Description:"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif"),
                     p("R implementation of the PhenoGraph algorithm A simple R implementation of the [PhenoGraph](http://www.cell.com/cell/abstract/S0092-8674(15)00637-6) algorithm, which is a clustering method designed for high-dimensional single-cell data analysis.
                       It works by creating a graph ('network') representing phenotypic similarities between cells by calclating the Jaccard coefficient between nearest-neighbor sets, and then identifying communities using the well known [Louvain method](https://sites.google.com/site/findcommunities/) in this graph."),
                     h4(strong("speed:"),icon("battery-quarter"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif")
                   ),
                   br(),
                   column(
                     width = 4,
                     sliderInput(inputId = "phenograph_K",label = "Knn K",min = 10,max = 1000,step = 10,value = 30)
                   ),
                   actionButton(inputId = "phenograph_run",label = "Calculate",style='margin-left:80%;width:15%;background:#31353e;color:#eeeeee')
                 ),
                 div(
                   id="Flowmeans_clust",
                   style="display:none",
                   h3(icon("desktop"),"FlowMeans"),
                   column(
                     width = 12,
                     h4(strong("Description:"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif"),
                     p("Finds a good fit to the data using k-means clustering algorithm. Then merges the adjacent dense spherical clusters to find non-spherical clusters."),
                     h4(strong("speed:"),icon("battery-three-quarters"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif")
                   ),
                   br(),
                   column(
                     width = 4,
                     selectInput(inputId ="flowMeans_Mahalanobis",label = "Mahalanobis",choices = c("TRUE","FALSE") ,selected = "TRUE"),
                     numericInput(inputId ="flowMeans_iter.max",label = "iter.max" ,value = 50,min = 1),
                     sliderInput(inputId = "flowMeans_K",label = "K value",min = 0,max = 150,step = 1,value = 0)
                   ),
                   column(
                     width = 4,
                     selectInput(inputId ="flowMeans_Standardize",label = "Standardize" ,choices = c("TRUE","FALSE") ,selected = "TRUE"),
                     selectInput(inputId ="flowMeans_Update",label = "Update" ,choices = c("Mahalanobis","Mean","None") ,selected = "Mahalanobis"),
                     sliderInput(inputId = "flowMeans_maxK",label = "maxK value",min = 0,max = 150,step = 1,value = 0)
                   ),
                   column(
                     width = 4,
                     selectInput(inputId ="flowMeans_addNoise",label = "addNoise",choices = c("TRUE","FALSE") ,selected = "TRUE"),
                     selectInput(inputId ="flowMeans_OrthagonalResiduals",label = "OrthagonalResiduals",choices = c("TRUE","FALSE") ,selected = "TRUE"),
                     sliderInput(inputId ="flowMeans_nstart",label = "nstart" ,value = 10,min = 1,max=200)
                   ),
                   actionButton(inputId = "flowMeans_run",label = "Calculate",style='margin-left:80%;width:15%;background:#31353e;color:#eeeeee')
                 ),
                 div(
                   id="SamSPECTRAL_clust",
                   style="display:none",
                   h3(icon("desktop"),"SamSPECTRAL"),
                   column(
                     width = 12,
                     h4(strong("Description:"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif"),
                     p("Given an FCS file as input, SamSPECTRAL first builds the communities to sample the data points.
                       Then, it builds a graph and after weighting the edges of the graph by conductance computation,
                       it is passed to a classic spectral clustering algorithm to find the spectral clusters.
                       The last stage of SamSPECTRAL is to combine the spectral clusters.
                       The resulting 'connected components' estimate biological cell populations in the data sample."),
                     h4(strong("speed:"),icon("battery-quarter"),sytle="font-size:18px,'Source Sans Pro','Helvetica Neue',Helvetica,Arial,sans-serif")
                   ),
                   br(),
                   column(
                     width = 4,
                     numericInput(inputId = "SamSPECTRAL_precision",label = "precision", min = 0 ,max =6,step=1,value = 6),
                     sliderInput(inputId = "SamSPECTRAL_k",label = "k.for_kmeans",min = 0,max = 300,step = 5,value = 0)
                   ),
                   column(
                     width = 4,
                     numericInput(inputId = "SamSPECTRAL_stabilizer",label = "stabilizer", min = 10 ,max =300,step=20,value = 100),
                     sliderInput(inputId = "SamSPECTRAL_separation",label = "separation.factor",min = 0.3,max = 2,step = 0.03,value = 0.39)
                   ),
                   column(
                     width = 4,
                     numericInput(inputId = "SamSPECTRAL_sigma",label = "normal.Sigma", min = 1 ,max =300,step=5,value = 200),
                     sliderInput(inputId = "SamSPECTRAL_number.of.clusters",label = "number.of.clusters",min = 0,max = 150,step = 5,value = 0)
                   ),
                   actionButton(inputId = "SamSPECTRAL_run",label = "Calculate",style='margin-left:80%;width:15%;background:#31353e;color:#eeeeee')
                 )
               )
           ),
           column(
             width = 2,
             mainPanel(
               width=12,
               style="background:#FFFFFF;overflow:auto;height:650px;border-radius:5px;padding:10px;border-style: solid;border-color:#aaa;border-width:1px",
               column(
                 width = 12,
                 h4(strong("Arguments:")),
                 div(
                   id="FlowSOM_desc",
                   style="display:none",
                   span(strong("max K:"),p("Maximum number of clusters to try out")),
                   br(),
                   span(strong("nClus:"),p(" Exact number of clusters to use. If not NULL, max will be ignored.")),
                   br(),
                   span(strong("method:"),p(" Clustering method to use, given as a string."))
                 ),
                 div(
                   id="DensityCut_desc",
                   style="display:none",
                   span(strong("Knn K:"),p("A integer to specify the number of neighbours in building the Knn graph.")),
                   br(),
                   span(strong("maxit:"),p("The maximum number of iteration allowed in density refinement, default to 50.")),
                   br(),
                   span(strong("eps:"),p(" The threshold in density refinement, default to 1e-5.")),
                   br(),
                   span(strong("alpha:"),p("The damping factor between 0 and 1, default to 0.90."))
                 ),
                 div(
                   id="Consboost_desc",
                   span(strong("K:"),p("Number of final clusters.")),
                   br(),
                   span(strong("sample_n:"),p("sampling the data for consensus clustering.")),
                   br(),
                   span(strong("n_Core:"),p("number of CPU threads be used for clustering"))
                 ),
                 div(
                   id="Flowmeans_desc",
                   style="display:none",
                   span(strong("K:"),p("Number of clusters. If set to '0' (default) the value will be estimated automatically.")),
                   br(),
                   span(strong("max K:"),p("Maximum number of clusters. If set to '0' (default) the value will be estimated automatically.")),
                   br(),
                   span(strong("Mahalanobis:"),p("Boolean value. If TRUE (default) mahalanobis distance will be used. Otherwised, euclidean distance will be used.")),
                   br(),
                   span(strong("Standardize:"),p("Boolean value. If TRUE (default) the data will be transformed to the [0,1] interval.")),
                   br(),
                   span(strong("addNoise:"),p("Boolean value. Determines if uniform noise must be added to the data to prevent singularity issues or not.")),
                   br(),
                   span(strong("iter.max:"),p("The maximum number of iterations allowed.")),
                   br(),
                   span(strong("Update:"),p("String value. If set to 'Mahalanobis' the distance function will be updated at each merging iteration with recalculating mahalanobis distances.")),
                   br(),
                   span(strong("OrthagonalResiduals:"),p("Boolean value, indicates if the residuals must be transformed to orthagonal distance or not.")),
                   br(),
                   span(strong("nstart:"),p("The number of random sets used for initialization."))
                 ),
                 div(
                   id = "Phenograph_desc",
                   style="display:none",
                   span(strong("K:"),p("integer; number of nearest neighbours (default:30).")),
                   br()
                 ),
                 div(
                   id = "SamSPECTRAL_desc",
                   style="display:none",
                   span(strong("precision:"),p("Determines the precision of computations. Setting it to 6 will work and increasing it does not improve results.")),
                   br(),
                   span(strong("stabilizer:"),p("The larger this integer is, the final results will be more stable because the underlying kmeans will restart many more times.")),
                   br(),
                   span(strong("normal.sigma:"),p("A scaling parameter that determines the 'resolution' in the spectral clustering stage. By increasing it, more spectral clusters are identified. This can be useful when 'small' population are aimed. See the user manual for a suggestion on how to set this parameter using the eigenvalue curve.")),
                   br(),
                   span(strong("k.for_kmeans:"),p("The number of clusters for running kmeans algorithm in spectral clustering. The default value of 'NA' leads to automatic estimation based on eigen values curve.")),
                   br(),
                   span(strong("separation.factor:"),p("This threshold controls to what extend clusters should be combined or kept separate.Normally, an appropriate value will fall in range 0.3-2.")),
                   br(),
                   span(strong("number.of.clusters:"),p("The default value is '0' which leads to computing the number of spectral clusters automatically, otherwise it can be a vector of integers each of which determines the number of spectral clusters. The output will contain a clustering resulting from each value..")),
                   br()
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
        h3("Cluster Label Generation",style="color:#111111"),
        style="background:#ffffff;padding:10px;border-radius:5px;",
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
             eChartOutput("Marker_line",height=680)
           ),
           actionButton(inputId = "back_clustering",label="back",width="20%",style="margin-top:1%;margin-left:10%"),
           actionButton(inputId = "skip1",label="Skip",width="20%",style="margin-top:1%;margin-left:10%"),
           actionButton(inputId = "label_button",label = "Label",width="20%",style="margin-top:1%;margin-left:10%;background:#31353e;color:#eeeeee")
         )
        )
      ),
      div(
        id="celltype_tree_div",
        h3("Clusters Table",style="color:#111111"),
        style="display:none;background:#ffffff;padding:10px;border-radius:5px;",
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
              visNetwork::visNetworkOutput("CellType_tree",height = 680,width = "100%")
            ),
            actionButton("Skip",label="Skip",width="25%",style="margin-top:1%;margin-left:17%"),
            actionButton("Next",label = "Next",width="25%",style="margin-top:1%;margin-left:15%;background:#31353e;color:#eeeeee")
          )
        )
      )
    ),
    tabItem(
     tabName ="ClusterinfoTable",
     div(style="background:#ffffff;padding:10px;border-radius:5px;",
      br(),
      fluidRow(
        column(width = 3),
        column(
          h4(strong("Clusters information:")),
          helpText("Revise the label for each cell population"),
          width = 6,
          uiOutput(outputId = "confirm_marker_choose"),
          plotOutput(outputId = "cell_type_confirm_plot",height = 330,width = "95%")

        )
      ),
      fluidRow(
        column(width = 3),
        column(
          width = 6,
          DT::dataTableOutput(outputId = "Cluster_info",width = "95%")
        )
      ),
      br(),
      actionButton("previous",label="Previous",width="15%",style="margin-bottom:1%;margin-left:30%"),
      actionButton("confirm",label="Confirm",width="15%",style="margin-bottom:1%;margin-left:10%;background:#31353e;color:#eeeeee")
     )
    ),

    #####=============== UI body Visualization ===============#####
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
             plotOutput(outputId = "scatter_plot",height = 670,width = "100%")
           ),
           column(
             width=3,
             box(
               width = 12,
               title = div(icon("cog"),"\t","Setting"),
               radioButtons(
                 inputId="red_choice",
                 label = "Dimensionality reduction methods:",
                 choices = c("t_SNE","LargeVis","PCA"),selected = "PCA",
                 inline = TRUE),
               br(),
               radioButtons(
                 inputId="Cluster_Marker",
                 label = "Choose a way to display:",
                 choices = c("Clusters","Markers"),selected = "Clusters",
                 inline = TRUE),
               br(),
               selectizeInput(
                 inputId = "scatter_colors",
                 label = "Color style:",
                 choices = c("set1","set2"),
                 selected = "set1"),
               br(),
               uiOutput("Marker_or_Cluster"),
               br()
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
             width=9,
             column(
               width = 6,
               plotOutput("MST_pie",height = 670)
             ),
             column(
               width = 5,
               plotOutput("MST_marker",height = 670)
             )
           ),
           column(
             width=3,
             box(
               width = 12,
               title = div(icon("cog"),"\t","Setting"),
               uiOutput("MST_panel",height=100),
               br()
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
            width = 9,
            plotOutput("Marker_heatmap",width = "100%",height = 670)
          ),
          column(
            width = 3,
            box(
              width = 12,
              title = div(icon("cog"),"\t","Setting"),
              column(
                width = 12,
                radioButtons(inputId = "heatmap_scale",label = "scale:",choices = c("row","column","none"),selected = "none",inline = TRUE),
                radioButtons(inputId = "heatmap_colorsets",label = "Color style:",choices = c("set1","set2"),selected = "set1",inline = TRUE),
                radioButtons(inputId = "heatmap_clustering_distance_rows",label = "clustering distance rows:",choices = c("correlation","euclidean"),selected = "euclidean",inline = TRUE),
                radioButtons(inputId = "heatmap_clustering_distance_cols",label = "clustering distance cols:",choices = c("correlation","euclidean"),selected = "euclidean",inline = TRUE),
                sliderInput(inputId = "heatmap_fontsize",label = "fontsize:",value = 10,min = 5,max = 15)
              )
            )
          )
         )
       ),

       ### DMT ###
       tabPanel(
         id="Marker_pop_vis",
         title="Population Marker",
         fluidRow(
          column(
            width = 6,
            eChartOutput("cell_type_marker_line",width = "100%",height = 670)
          ),
          column(
            width = 6,
            visNetwork::visNetworkOutput("DMT",width = "100%",height = 670)
          )
         )
       ),
       tabPanel(
         id = "result_report",
         title = "Final report",
         column(
           width = 3
         ),
         column(
           width = 6,
           div(
             radioButtons(inputId = "choose_report",label = "select report",choices = c("Project","Cluster"),selected = "Project",inline = TRUE),
             div(
               id="Project_report_div",
               h4(strong("Project Report:")),
               DT::dataTableOutput(outputId = "Project_report",width = "100%")
             ),
             div(
               id="Cluster_report_div",
               style="display:none",
               h4(strong("Cluster Report")),
               helpText("If there is no autmatic labeling method be used, This table will not show up."),
               DT::dataTableOutput(outputId = "Cluster_report",width = "100%")
             )
           )
         )
       )
     ),
     absolutePanel(
       width = "230px",
       draggable = TRUE,
       bottom = 30,
       right = 30,
       div(
        style="width:100%;background:#eeeeee;border-radius:5px;",
        h4(icon("cog"),"\t","Setting Panel",style="background:#121212;color:#efefef;padding:5px"),
        div(
          style="width:100%;background:#eeeeee;padding:15px;text-align: center",
          radioButtons(inputId = "IDexchangeLabel",label = "ID transform into Cell Label:",choices = c("ID","Cell_Label"),inline = TRUE),
          actionButton(inputId = "back_labeling",label = "Relabel",width="70%"),
          downloadButton(outputId = "DownloadRdata",
                         label = "Download Rdata",
                         style="width:70%;margin-top:5%;background:#333333;color:#eeeeee;border-radius:5px")
        )
      )
     )
    )
   )
  )
 ) #dashbody end
)  #UI end



###===========================shiny servers part================================####

server = function(input,output,session){
 options(shiny.maxRequestSize=250*1024^2)
 options(rgl.useNULL = TRUE)

 ###=========================== step one ===========================####

 # vector for saving choosen markers #
 dy_markers <- reactiveValues(data=NULL)

 # vector for refreshing the graph
 # we will use refresh$dynamic+=1 and refresh$dynamic to refresh pictures.
 refresh <- reactiveValues(dynamic=1)

 # a vector for saving file #
 file <- reactiveValues(data=NULL)

 # file loading setting #
 observeEvent(input$fcs,{
   # saving result as cytoSet
   inFile <- input$fcs
   file$data<-read.FCS(inFile$datapath)
   cyto<<-initiflow(fcs.data = exprs(file$data),projectname = input$Project_name1)
   cyto@filename<<-inFile$name
   cyto@version<<-packageVersion("cytosee")
   cyto@autoLabel<<-FALSE
   cyto@preprocess<<-"None"
   showElement("Step1_body",animType ="slide",anim = TRUE,time = 0.5)
   showElement("Step1_bar",animType ="slide",anim = TRUE,time = 0.5)
   hideElement(selector = ".Introduction")
   hideElement(selector = ".InputFile")
   showElement(selector = ".Preview")

 })

 observeEvent(input$Rdata,{
   inFile <- input$Rdata
   load(inFile$datapath)
   if(is.null(cyto@ClusterID)){
     shinyalert("Invaild input file !",
                type = "error",
                confirmButtonCol = "#3c8dbc")
      return()
   }
   else{
     hideElement(selector = ".Introduction")
     hideElement(selector = ".InputFile")
     showElement(selector=".Step3")
     }
 })


 ## plot bin set
 output$Plot_bin<-renderUI({
   if(is.null(file$data)){
     return()
   }
   fluidRow(
     column(
      width=6,
      sliderInput(inputId = "bin_val",
                  label = "bin value",
                  min=0,
                  max=512,
                  value = 256,
                  step=16)
     ),
     column(
       width=6,
       sliderInput(inputId = "scale_val",
                   label = "scale value",
                   min=10,
                   max=500,
                   value = 30,
                   step=10)
     )
   )
 })

 #####data Preview hex
 output$hex_plot<-renderPlot({
   data<- file$data
   if(is.null(input$transform)){
     return()
   }

   if (input$transform=="log10"){
     ### check data ###
     # make sure the transform methods can be used
     X=data@exprs
     for(i in 1:length(X[1,])){
       minimal=min(X[,i])
       if(minimal<=0 && tolower(colnames(X)[i])!="time"){
         X[,i]=X[,i]+abs(minimal)+0.1
       }
     }
     data@exprs=X
     tran <- logTransform(transformationId="log10-transformation", logbase=10, r=1, d=1)
     channeltouse=input$choose_marker_trans
     tranlist_log <- transformList(channeltouse, tran)
     data<-transform(data, tranlist_log)
     cyto@transform_method<<-"log10"
   }

   else if (input$transform=="arcsinh"){
     tran <- arcsinhTransform("arcsinh-transformation",a = 1,b=1,c=0)
     channeltouse=input$choose_marker_trans
     tranlist_arcsinh <- transformList(channeltouse, tran)
     data<-transform(data, tranlist_arcsinh)
     cyto@transform_method<<-"arcsinh"
   }

   else{
     cyto@transform_method<<-"None"
   }
   cyto@transform_paras<<-input$choose_marker_trans
   data_trans$data<-data
   scale=input$scale_val
   cyto@fcs.data<<-as.data.frame(exprs(data))
   cyto@event.use<<-c(1:length(cyto@fcs.data))
   data<-cyto@fcs.data
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
     scale_fill_gradientn(colours = c("black","darkred","red","orange","yellow","white"),
                          limits=c(0,scale),
                          name="Count of cells\n")
 })

 ## Data Preview density plot  ##
 output$X_density<-renderPlot({
   if(is.null(data_trans$data)){
     return()
   }
   data=as.data.frame(exprs(data_trans$data))
   X=data[input$plotx]
   diff=abs(range(X)[2]-range(X)[1])
   bin=diff/200
   ggplot(NULL,aes(x=X))+
     geom_histogram(binwidth =bin,fill="#c7b3e5",colour="#c7b3e5")+
     ggtitle("Marker Density Map\n")+
     theme_bw()+
     theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
     theme(panel.border  = element_blank())+
     theme(axis.line = element_line(colour = "black"))+xlab(input$plotx)
 })

 output$Y_density<-renderPlot({
   if(is.null(data_trans$data)){
     return()
   }
   data=as.data.frame(exprs(data_trans$data))
   Y=data[input$ploty]
   diff=abs(range(Y)[2]-range(Y)[1])
   bin=diff/200
   ggplot(NULL,aes(x=Y))+
     geom_histogram(binwidth =bin,fill="#daf9ca",colour="#daf9ca")+
     ggtitle("Marker Density Map\n")+
     theme_bw()+
     theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
     theme(panel.border  = element_blank())+
     theme(axis.line = element_line(colour = "black"))+xlab(input$ploty)
 })



 ###### file description #
 output$file <- DT::renderDataTable({
   inFile <- file$data
   if (is.null(inFile))
     return(NULL)
   fcs=file$data
   # define the base view table for files
   table<-data.frame(as.data.frame(fcs@parameters@data$name),
                     as.data.frame(fcs@parameters@data$desc),
                     row.names = NULL)

   colnames(table)<-c("name","desc")

   DT::datatable(
     table,
     class = "hover",
     style = "bootstrap",
     options = list(dom="rpit",paging=FALSE,scrollY=250),
     selection = list(mode="multiple",
                      selected=r_all_choice$data,
                      target="row")
   )
 })

 # choose the marker for data transformation
 output$Choose_transform_channel <- renderUI({
   inFile <- file$data
   if (is.null(inFile))
     return(NULL)
   fcs=file$data
   name<-as.character(fcs@parameters@data$name)
   checkboxGroupInput(inputId = "choose_marker_trans",label = "",choices = name,selected = name,inline = TRUE)
 })

 # choose the display markers
 output$XY_val <- renderUI({
   fcs<-file$data
   if (is.null(fcs)){
     return()
   }
   tagList(
     fluidRow(
       column(
         width = 6,
         selectInput("plotx", "plot-x",colnames(fcs),selected = colnames(fcs)[1])
         ),
       column(
         width = 6,
         selectInput("ploty", "plot-y",colnames(fcs),selected = colnames(fcs)[2]))
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

 # function of button "back_homepage"
 observeEvent(input$back_homepage,{
   shinyalert(
     title = "GO back?",
     text = "Click 'Confirm' will go back to the homepage and import data will all be deleted from cytosee object",
     type = "warning",
     confirmButtonCol = "#3c8dbc",
     confirmButtonText = "Confirm",
     showCancelButton = TRUE,
     callbackR = function(x) { if(x != FALSE){
       cyto<<-NA
       hideElement("Step1_body",animType ="slide",anim = TRUE,time = 0.3)
       hideElement("Step1_bar",animType ="slide",anim = TRUE,time = 0.3)
       showElement(selector = ".Introduction")
       showElement(selector = ".InputFile")
       hideElement(selector = ".Preview")} }
   )
 })

 # function of button "all" #
 r_all_choice=reactiveValues(data=NULL)
 observeEvent(input$all_choice,{
   refresh$dynamic
   r_all_choice$data=NULL
   r_all_choice$data=c(1:length(file$data@parameters@data$name))
 })

 # function of button "Next"#
 observeEvent(input$previewer_next,{
   inFile<-input$fcs
   refresh$dynamic=refresh$dynamic+0.01
   if (is.null(inFile)){
     return("A .fcs is needed")
   }
   fcs=read.FCS(inFile$datapath)
   channel_id = input$file_rows_selected
   dy_markers$data=colnames(exprs(fcs))[as.integer(channel_id)]
   if(length(dy_markers$data)==0){
     shinyalert("No Channel was chosen, at least one Channel is required!",
                type = "error",
                confirmButtonCol = "#3c8dbc")
     return()
   }
   hideElement(selector = ".Step1")
   hideElement(id="Step1_bar")
   showElement(selector = ".Preprocess")
   hideElement(id="flowAI_UI")
   showElement(id="choose_QC")
   cyto@channel.use<<-channel_id
 })


#####================ preprocess step ================#####

# function of skip
observeEvent(input$select_QC_skip,{
   refresh$dynamic=refresh$dynamic+1
   data=data_trans$data
   if(is.null(data)){
     return()
   }
   hideElement(selector = ".Preprocess",animType ="fade",anim = TRUE,time = 0.3)
   showElement(selector = ".Step2",animType ="slide",anim = TRUE,time = 1)
   showElement("Step2_bar")
   cyto@event.use<<-c(1:length(cyto@fcs.data[,1]))
   cyto@preprocess<<-"NaN"
})

# function of next
observeEvent(input$select_QC_next,{
  refresh$dynamic=refresh$dynamic+1
  hideElement(id= "choose_QC")
  if(input$select_QC=="flowAI-manual"){
    showElement(id="flowAI_UI")
  }
  else if(input$select_QC=="flowAI-auto"){
    withProgress(message = 'Run QC method',
                 detail = 'This may take a while...',
                 value = 0, {
                   data<-data_trans$data
                   incProgress(5/10,message="Running flowAI method")
                   re<-cytosee_flowAI(data,html_report = FALSE,folder_results = FALSE,output = 3,fcs_QC = FALSE,mini_report = FALSE)
                   cyto@event.use<<-setdiff(1:length(cyto@fcs.data[,1]),re$`0`)
                   cyto@preprocess<<-"flowAI-manual"
                   incProgress(10/10,message="Completed!")
                 })
    hideElement(selector = ".Preprocess",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step2",animType ="slide",anim = TRUE,time = 1)
    showElement("Step2_bar")
  }
})

# PreProcess button for go back to preview page #
 observeEvent(input$back_preview,{
   shinyalert(
     title = "GO back?",
     text = "Click 'Confirm' will go back to the preview page and data produced from preview will all be deleted from cytosee object",
     type = "warning",
     confirmButtonCol = "#3c8dbc",
     confirmButtonText = "Confirm",
     showCancelButton = TRUE,
     callbackR = function(x) { if(x != FALSE){
       hideElement(selector = ".Preprocess",animType ="slide",anim = TRUE,time = 0.3)
       showElement(selector = ".Step1",animType ="slide",anim = TRUE,time = 0.3)
       showElement("Step1_bar")
       cyto@channel.use<<-numeric()
       cyto@transform_method<<-""
       }
      }
   )
 })
 ## flowAI-manual Servers ##

 # flowAI button for Skipping #
 observeEvent(input$AI_Skip,{
   refresh$dynamic=refresh$dynamic+1
   data=data_trans$data
   if(is.null(data)){
     return()
   }
   hideElement(selector = ".Preprocess",animType ="fade",anim = TRUE,time = 0.3)
   showElement(selector = ".Step2",animType ="slide",anim = TRUE,time = 1)
   showElement("Step2_bar")
   cyto@event.use<<-c(1:length(cyto@fcs.data[,1]))
   cyto@preprocess<<-"NaN"
 })


  # flowAI button for confirming #
 observeEvent(input$AI_Confirm,{
   refresh$dynamic=refresh$dynamic+1
   data=data_trans$data
   if(is.null(data)){
     return()
   }
   hideElement(selector = ".Preprocess",animType ="fade",anim = TRUE,time = 0.3)
   showElement(selector = ".Step2",animType ="slide",anim = TRUE,time = 1)
   showElement("Step2_bar")
   re<-checkRes()[[9]]
   cyto@event.use<<-re
   cyto@preprocess<<-"flowAI-manual"
 })

  ## load flowset data
  set <- reactive({
    refresh$dynamic
    isolate({fcsFiles <- input$fcs
    if (is.null(fcsFiles))
      return(NULL)
    set <- data_trans$data
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
    if(is.null(set()) || is.null(cellCheck())){
      return(NULL)
    }
    flowSignalData <- cellCheck()[[2]]
    fsp <- flow_signal_plot(flowSignalData, input$signalBinSlider[1], input$signalBinSlider[2])
    print(fsp)
  })

  output$flowSignalPlotUI <- renderUI({
    if(is.null(set()) || is.null(cellCheck())){
      return(NULL)
    }
    plotOutput("flowSignalPlot",height = 45*length(colnames(cyto@fcs.data)))
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
                flowRatePerc, flowSignalPerc, flowMarginPerc, QCfcs,goodCellIDs))
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

    })

###================================Step 2 Cell clustering ================================####

  ##### User submit clusters #####
  observeEvent(input$import_cluster,{
    refresh$dynamic=refresh$dynamic+0.01
    inFile <- input$import_cluster
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   file=flowFrame(exprs = as.matrix(cyto@fcs.data[cyto@event.use,]))
                   incProgress(1/10,message="Import ClusterID ...")
                   ClusterID = read.csv(file=inFile$datapath,header = TRUE)
                   ClusterID = as.data.frame(ClusterID[cyto@event.use,])
                   data=cbind(ClusterID,cyto@fcs.data[cyto@event.use,][,cyto@channel.use])
                   cyto@label<<-ClusterID
                   Cluster_id <- unique(cyto@label)
                   cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                   cyto@ClusterID<<-ClusterID
                   cyto@clust_method<<-"External"
                   incProgress(1/10,message="Build DMT...")
                   cyto@dmt<<-cytosee_DMT(data)
                   incProgress(1/10,message="Build MST...")
                   if(input$pre_mst_method=="flowSOM"){
                     file=flowFrame(exprs = as.matrix(cyto@fcs.data[cyto@event.use,]))
                     fsom=BuildSOM(ReadInput(file),colsToUse = cyto@channel.use)
                     MST=BuildMST(fsom)
                   }
                   cyto@mst<<-MST
                   incProgress(7/10,message="Reduce the Dimension...")
                   cyto@dim.red<<-cytosee_DR(cyto@fcs.data[cyto@event.use,][,cyto@channel.use],method = input$pre_reduce_method)
                 })
    shinyalert(
      title = "Complete!",
      type = "success",
      confirmButtonCol = "#3c8dbc",
      confirmButtonText = "OK"
    )
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })


  ## cluster choosing event ##
  observeEvent(input$select_clust_method,{
    clust_name=c("Consboost","DensityCut","Flowmeans","FlowSOM","Phenograph","SamSPECTRAL")
    for(i in 1:6){
        hideElement(sprintf("%s_clust",clust_name[i]))
        hideElement(sprintf("%s_desc",clust_name[i]))
    }
    showElement(sprintf("%s_clust",input$select_clust_method))
    showElement(sprintf("%s_desc",input$select_clust_method))
  })

  ## button for going back to preprocess page ##
  observeEvent(input$back_preprocess,{
    shinyalert(
      title = "GO back?",
      text = "Click 'Confirm' will go back to the preview page and data produced from preprocess will all be deleted from cytosee object",
      type = "warning",
      confirmButtonCol = "#3c8dbc",
      confirmButtonText = "Confirm",
      showCancelButton = TRUE,
      callbackR = function(x) { if(x != FALSE){
        showElement(selector = ".Preprocess",animType ="slide",anim = TRUE,time = 0.3)
        hideElement(selector = ".Step2")
        hideElement(id="flowAI_UI")
        hideElement("Step2_bar")
        cyto@event.use<<-numeric()
        cyto@preprocess<<-""
        }
      }
    )
  })

  ##### flowSOM #####
  observeEvent(input$flowSOM_run,{
    refresh$dynamic=refresh$dynamic+0.01
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                  file=flowFrame(exprs = as.matrix(cyto@fcs.data[cyto@event.use,]))
                  incProgress(1/10,message="Run flowSOM...")
                  if(input$flowSOM_nclus=="NULL"){
                    nclus=NULL
                  }
                  else{
                    nclus=as.numeric(input$flowSOM_nclus)
                  }
                  cyto <<- cytosee_flowSOM(cyto,K = input$flowSOM_K,nClus=nclus,method=input$flowSOM_metaclust_method)
                  data=cbind(cyto@ClusterID,cyto@fcs.data[cyto@event.use,][,cyto@channel.use])
                  cyto@label<<-cyto@ClusterID
                  Cluster_id <- unique(cyto@label)
                  cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                  cyto@clust_method<<-"flowSOM"
                  incProgress(1/10,message="Build DMT...")
                  cyto@dmt<<-cytosee_DMT(data)
                  incProgress(1/10,message="Build MST...")
                  if(input$pre_mst_method=="flowSOM"){
                    file=flowFrame(exprs = as.matrix(cyto@fcs.data[cyto@event.use,]))
                    fsom=BuildSOM(ReadInput(file),colsToUse = cyto@channel.use)
                    MST=BuildMST(fsom)
                  }
                  cyto@mst<<-MST
                  incProgress(6/10,message="Reduce the Dimension...")
                  cyto@dim.red<<-cytosee_DR(cyto@fcs.data[cyto@event.use,][,cyto@channel.use],method = input$pre_reduce_method)
                  incProgress(1/10,message="Complete!")
                 })
    shinyalert(
      title = "Complete!",
      type = "success",
      confirmButtonCol = "#3c8dbc",
      confirmButtonText = "OK"
    )
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })

  ##### DensityCut #####
  observeEvent(input$DensityCut_run,{
    refresh$dynamic=refresh$dynamic+0.01
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   file= cyto@fcs.data[cyto@event.use,][,cyto@channel.use]
                   incProgress(1/10,message="Run DesnityCut")
                   cyto <<- cytosee_DensityCut(cyto,K = input$DensityCut_K,alpha=input$DensityCut_alpha,maxit=input$DensityCut_maxit,eps=input$DensityCut_eps)
                   data=cbind(cyto@ClusterID,cyto@fcs.data[cyto@event.use,][,cyto@channel.use])
                   cyto@label<<-cyto@ClusterID
                   Cluster_id <- unique(cyto@label)
                   cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                   cyto@clust_method<<-"DensityCut"
                   incProgress(1/10,message="Build DMT...")
                   cyto@dmt<<-cytosee_DMT(data)
                   incProgress(1/10,message="Build MST...")
                   if(input$pre_mst_method=="flowSOM"){
                     file=flowFrame(exprs = as.matrix(cyto@fcs.data[cyto@event.use,]))
                     fsom=BuildSOM(ReadInput(file),colsToUse = cyto@channel.use)
                     MST=BuildMST(fsom)
                   }
                   cyto@mst<<-MST
                   incProgress(6/10,message="Reduce the Dimension...")
                   cyto@dim.red<<-cytosee_DR(cyto@fcs.data[cyto@event.use,][,cyto@channel.use],method = input$pre_reduce_method)
                   incProgress(1/10,message="Complete!")
                 })
    shinyalert(
      title = "Complete!",
      type = "success",
      confirmButtonCol = "#3c8dbc",
      confirmButtonText = "OK"
    )
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })


  ##### Consboost #####
  observeEvent(input$Consboost_run,{
    refresh$dynamic=refresh$dynamic+0.01
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   incProgress(1/10,message="Run Consboost")
                   cyto<<-cytosee_consboost(cyto,K=input$Consboost_K,sample_n=input$Consboost_sample_n,n_core=input$Consboost_n_core)
                   data=cbind(cyto@ClusterID,cyto@fcs.data[cyto@event.use,][,cyto@channel.use])
                   cyto@label<<-cyto@ClusterID
                   Cluster_id <- unique(cyto@label)
                   cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                   cyto@clust_method<<-"Consboost"
                   incProgress(1/10,message="Build DMT...")
                   cyto@dmt<<-cytosee_DMT(data)
                   incProgress(1/10,message="Build MST...")
				   file= cyto@fcs.data[cyto@event.use,][,cyto@channel.use]
                   if(input$pre_mst_method=="flowSOM"){
                     file=flowFrame(exprs = as.matrix(cyto@fcs.data[cyto@event.use,]))
                     fsom=BuildSOM(ReadInput(file),colsToUse = cyto@channel.use)
                     MST=BuildMST(fsom)
                   }
                   cyto@mst<<-MST
                   incProgress(6/10,message="Reduce the Dimension...")
                   cyto@dim.red<<-cytosee_DR(cyto@fcs.data[cyto@event.use,][,cyto@channel.use],method = input$pre_reduce_method)
                   incProgress(1/10,message="Complete!")
                 })
    shinyalert(
      title = "Complete!",
      type = "success",
      confirmButtonCol = "#3c8dbc",
      confirmButtonText = "OK"
    )
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })


  ##### FlowMeans #####
  observeEvent(input$flowMeans_run,{
    refresh$dynamic=refresh$dynamic+0.01
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   file= cyto@fcs.data[cyto@event.use,][,cyto@channel.use]
                   incProgress(1/10,message="Run FlowMeans")
        				   if(input$flowMeans_maxK==0){
        					    maxK=NA
        				   }
        				   else{
        					    maxK=input$flowMeans_maxK
        				   }
        				   if(input$flowMeans_K==0){
            					K=NA
        				   }
        				   else{
        				      K=input$flowMeans_K
        				   }
                   cyto <<- cytosee_flowMeans(
                     cyto,
                     MaxN=maxK,
                     NumC=K,
                     Mahalanobis=as.logical(input$flowMeans_Mahalanobis),
                     Standardize=as.logical(input$flowMeans_Standardize),
                     addNoise=as.logical(input$flowMeans_addNoise),
                     iter.max=input$flowMeans_iter.max,
                     Update=input$flowMeans_Update,
                     OrthagonalResiduals=as.logical(input$flowMeans_OrthagonalResiduals),
                     nstart=input$flowMeans_nstart
				          )
                   data=cbind(cyto@ClusterID,cyto@fcs.data[cyto@event.use,][,cyto@channel.use])
                   cyto@label <<- cyto@ClusterID
                   Cluster_id <- unique(cyto@label)
                   cyto@ID2CL <<- data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                   cyto@clust_method <<- "FlowMeans"
                   incProgress(1/10,message="Build DMT...")
                   cyto@dmt <<- cytosee_DMT(data)
                   incProgress(1/10,message="Build MST...")
                   if(input$pre_mst_method=="flowSOM"){
                     file = flowFrame(exprs = as.matrix(cyto@fcs.data[cyto@event.use,]))
                     fsom = BuildSOM(ReadInput(file),colsToUse = cyto@channel.use)
                     MST = BuildMST(fsom)
                   }
                   cyto@mst<<-MST
                   incProgress(6/10,message="Reduce the Dimension...")
                   cyto@dim.red<<-cytosee_DR(cyto@fcs.data[cyto@event.use,][,cyto@channel.use],method = input$pre_reduce_method)
                   incProgress(1/10,message="Complete!")
                 })
    shinyalert(
      title = "Complete!",
      type = "success",
      confirmButtonCol = "#3c8dbc",
      confirmButtonText = "OK"
    )
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })


  ##### Phenograph #####
  observeEvent(input$phenograph_run,{
    refresh$dynamic=refresh$dynamic+0.01
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   file= cyto@fcs.data[cyto@event.use,][,cyto@channel.use]
                   incProgress(1/10,message="Run Phenograph")
                   cyto <<- cytosee_phenograph(cyto,K=input$phenograph_K)
                   data=cbind(cyto@ClusterID,cyto@fcs.data[cyto@event.use,][,cyto@channel.use])
                   cyto@label<<-cyto@ClusterID
                   Cluster_id <- unique(cyto@label)
                   cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                   cyto@clust_method<<-"Phenograph"
                   incProgress(1/10,message="Build DMT...")
                   cyto@dmt<<-cytosee_DMT(data)
                   incProgress(1/10,message="Build MST...")
                   if(input$pre_mst_method=="flowSOM"){
                     file=flowFrame(exprs = as.matrix(cyto@fcs.data[cyto@event.use,]))
                     fsom=BuildSOM(ReadInput(file),colsToUse = cyto@channel.use)
                     MST=BuildMST(fsom)
                   }
                   cyto@mst<<-MST
                   incProgress(6/10,message="Reduce the Dimension...")
                   cyto@dim.red<<-cytosee_DR(cyto@fcs.data[cyto@event.use,][,cyto@channel.use],method = input$pre_reduce_method)
                   incProgress(1/10,message="Complete!")
                 })
    shinyalert(
      title = "Complete!",
      type = "success",
      confirmButtonCol = "#3c8dbc",
      confirmButtonText = "OK"
    )
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })


  ##### SamSPECTRAL #####
  observeEvent(input$SamSPECTRAL_run,{
    refresh$dynamic=refresh$dynamic+0.01
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   incProgress(1/10,message="Run SamSPECTRAL...")
        				   if(input$SamSPECTRAL_k==0){
        					 k = "NA"
        				   }
        				   else{
        					 k = input$SamSPECTRAL_k
        				   }
        				   if(input$SamSPECTRAL_number.of.clusters==0){
        					 number.of.clusters="NA"
        				   }
        				   else{
        					 number.of.clusters=input$SamSPECTRAL_number.of.clusters
        				   }
                   cyto<<-cytosee_SamSPECTRAL(cyto,
                      											  sigma=input$SamSPECTRAL_sigma,
                                              separation.factor=input$SamSPECTRAL_separation,
                      											  precision=input$SamSPECTRAL_precision,
                      											  k.for_kmeans=k,
                      											  stabilizer=input$SamSPECTRAL_stabilizer,
                      											  number.of.clusters=number.of.clusters
					         )
                   cyto@ClusterID[which(is.na(cyto@ClusterID)),]="Unknow"
                   cyto@ClusterID<<-cyto@ClusterID
                   data=cbind(cyto@ClusterID,cyto@fcs.data[cyto@event.use,][,cyto@channel.use])
                   cyto@label<<-cyto@ClusterID
                   Cluster_id <- unique(cyto@label)
                   cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
                   cyto@clust_method<<-"SamSPECTRAL"
                   incProgress(1/10,message="Build DMT...")
                   cyto@dmt<<-cytosee_DMT(data)
                   incProgress(1/10,message="Build MST...")
                   if(input$pre_mst_method=="flowSOM"){
                     file=flowFrame(exprs = as.matrix(cyto@fcs.data[cyto@event.use,]))
                     fsom=BuildSOM(ReadInput(file),colsToUse = cyto@channel.use)
                     MST=BuildMST(fsom)
                   }
                   cyto@mst<<-MST
                   incProgress(6/10,message="Reduce the Dimension...")
                   cyto@dim.red<<-cytosee_DR(cyto@fcs.data[cyto@event.use,][,cyto@channel.use],method = input$pre_reduce_method)
                   incProgress(1/10,message="Complete!")
                 })
    shinyalert(
      title = "Complete!",
      type = "success",
      confirmButtonCol = "#3c8dbc",
      confirmButtonText = "OK"
    )
    hideElement(selector = ".Step2",animType ="fade",anim = TRUE,time = 0.3)
    showElement(selector = ".Step3",animType ="slide",anim = TRUE,time = 1)
  })

  output$File_info<-renderUI({
    refresh$dynamic
    mainPanel(
      width=12,
      style="background:#FFFFFF;height:650px;overflow:auto;max-height:800px;border-radius:5px;padding:10px;border-style: solid;border-color:#aaa;border-width:1px",
      h4(strong("File information:")),
      column(
        width = 12,
        h5(strong("file name:")),
        span(cyto@filename),
        br(),
        h5(strong("cell num:")),
        span(length(cyto@event.use)),
        br(),
        h5(strong("low quality cell:")),
        span(length(cyto@fcs.data[,1])-length(cyto@event.use)),
        br(),
        h5(strong("transformation method:")),
        span(cyto@transform_method),
        br(),
        h5(strong("clustering channels:")),
        div(
          style="background : #eee",
          dataTableOutput("clust_channel_view")
        )
      )
     )
  })

  output$clust_channel_view<-renderDataTable({
    refresh$dynamic
    datatable=data.frame(colnames(cyto@fcs.data[cyto@channel.use]))
    colnames(datatable)="Channels"
    datatable
  },
    option=list(scrollY=250,dom="rpt",paging = FALSE,ordering=FALSE,searching=FALSE)
  )

#####================================Step 3 Labeling and Visualization================================#####


  ###=================== Label the clusters ===================###
  output$alias<-renderUI({
    refresh$dynamic
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
                    value = markers[i],
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
      height:500px',
      taglist
    )
  })

  ## button for going back to clustering step ##
  observeEvent(input$back_clustering,{
    shinyalert(
      title = "GO back?",
      text = "Click 'Confirm' will go back to the clustering step and data produced from clustering will all be deleted from cytosee object",
      type = "warning",
      confirmButtonCol = "#3c8dbc",
      confirmButtonText = "Confirm",
      showCancelButton = TRUE,
      callbackR = function(x) { if(x != FALSE){
        cyto@dmt <<- list()
        cyto@mst <<- list()
        cyto@label <<- list()
        cyto@ClusterID <<- list()
        cyto@dim.red <<- list()
        showElement(selector = ".Step2",animType ="slide",anim = TRUE,time = 0.3)
        hideElement(selector = ".Step3",animType ="fade",anim = TRUE,time = 0.3)
        }
      }
    )
  })


  ## button for Skipping label step ##
  r_marker_value=reactiveValues(data=c())
  observeEvent(input$skip1,{
    refresh$dynamic=refresh$dynamic+0.01
    Cluster_id <- unique(cyto@label)
    cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
    updateTabItems(session,"result_display",selected="ClusterinfoTable")
  })

  ## label button ##
  observeEvent(input$label_button,{
    refresh$dynamic=refresh$dynamic+0.01
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
      if(eval(parse(text = sprintf("input$%s !=''",marker_id)))){
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
      shinyalert("Input area shouldn't be empty!",
                 type = "error",
                 confirmButtonCol="#3c8dbc")
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
      word_p=paste0(word,">")
      word_m=paste0(word,"<")
      for(i in c(1:length(re))){
        store=re[i]
        # first and second symbol will determine this marker is "+" or "-"
        first=""
        second=""
        cnt=0
        while(length(grep(pattern = word_m,x = store,fixed = TRUE))!=0 || length(grep(pattern =word_p,x = store,fixed = TRUE))!=0 ){
          if(length(grep(word_p,store))!=0){
            cnt=cnt+1
            store=sub(pattern=word_p,replacement="",store,perl = FALSE)
            if(cnt==1){
              first="plus"
            }
            else if(cnt==2){
              second="plus"
            }
          }
          if(length(grep(word_m,store))!=0){
            cnt=cnt+1
            store=sub(pattern = word_m,replacement="",store,perl = FALSE)
            if(cnt==1){
              first="minus"
            }
            else if(cnt==2){
              second="minus"
            }
          }
        }
        if(cnt==0){
          next()
        }
        if(first=="plus"&&second=="plus"){
          re[i]=paste0(store,word,"++")
        }
        else if(first=="plus"){
          re[i]=paste0(store,word,"+")
        }
        else if(first=="minus"&&second=="minus"){
          re[i]=paste0(store,word,"--")
        }
        else if(first=="minus"){
          re[i]=paste0(store,word,"-")
        }
      }
    }
    CL_label=list()
    tryCatch({
      withProgress(
        message = "Run cytosee_LocCL",
        value=0,{
          for(i in 1:length(re)){
            incProgress(1/length(re))
            result=cytosee_LocCL(MarkerList = re[i],MaxHitsPht = input$CL_MaxHitsPht)
            ClusterID=Cluster_m2c[i]
            CL_label[[i]]=list("result"=result,"ClusterID"=ClusterID,"PhenoType"=re[i])
          }
        })
      cyto@CL_label<<-CL_label
      cyto@autoLabel<<-TRUE
      hideElement(id="change_marker_name")
      showElement(id="celltype_tree_div")
      },
      error=function(e){
        shinyalert("Something wrong with your query markers!",
                   type = "error",
                   confirmButtonCol="#3c8dbc")
        CL_label<-list()
        return()
    }
    )
  })



  ## Marker_line ##
  output$Marker_line<-renderEChart({
    refresh$dynamic
    cl_average_vec<-getAvg(fcsData=cyto@fcs.data[cyto@event.use,][,cyto@channel.use],ClusterID=cyto@label)
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
      setLegend(pos=6) %>%
      setTheme('macarons') %>%
      setToolbox(language = "en") %>%
      setTitle(title = "Marker expression level in each clusters",pos = 12)
  })


  # helper function for making checkbox in CL label datatable
  shinyInput_CL = function(FUN, len, id, ...){
    refresh$dynamic
    inputs = character(len)
    for (i in seq_len(len)) {
      inputs[i] = as.character(FUN(paste0(id, i),paste0(cyto@CL_label[[i]]$result$Labels),...))
    }
    inputs
  }

  # helper function for making checkbox in Cluster report datatable
  shinyInput_Cluster = function(FUN, len, id,Label, ...){
    refresh$dynamic
    inputs = character(len)
    for (i in seq_len(len)) {
      inputs[i] = as.character(FUN(paste0(id, i),paste0(Label[i]),...))
    }
    inputs
  }


  # helper function for reading checkbox in datatable
  shinyValue = function(id, len) {
    refresh$dynamic
    unlist(lapply(seq_len(len), function(i) {
      value = input[[paste0(id, i)]]
      if (is.null(value))
        NA
      else
        value
    }))
  }

  #####=================== clusters label page ===================#####
  output$Clusters_table<-DT::renderDataTable({
    refresh$dynamic
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
    refresh$dynamic
    id=input$Clusters_table_rows_selected
    if(is.null(id)){
      return("you can choose a cluster to display")
    }
    Nodes<-cyto@CL_label[[id]]$result$Nodes
    Links<-cyto@CL_label[[id]]$result$Links
    visNetwork(Nodes, Links) %>%
      visEdges(arrows = "from",color="#cccccc") %>%
      visNodes(size=18) %>%
      visExport(type="png",name="export-tree",label="Export as png") %>%
      visOptions(highlightNearest = list(enabled=TRUE,algorithm="hierarchical",
                                         degree=list(from=0,to=50))) %>%
      visHierarchicalLayout(nodeSpacing=150,
                            levelSeparation = 110,direction = "DU",
                            blockShifting = FALSE,sortMethod = "directed")
  })

  ##### Skip #####
  observeEvent(input$Skip,{
    refresh$dynamic=refresh$dynamic+0.01
    Cluster_id <- unique(cyto@label)
    cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
    updateTabItems(session,"result_display",selected="ClusterinfoTable")
  })

  ##### Next button #####
  observeEvent(input$Next,{
    refresh$dynamic=refresh$dynamic+0.01
    Cluster_id <- unique(cyto@label)
    cyto@ID2CL <<- data.frame(ID=Cluster_id,CD_Label = shinyValue("selecter_", nrow(Cluster_id)))
    updateTabItems(session,"result_display",selected="ClusterinfoTable")
  })


  #####=================== Cluster_info page ===================#####
  output$Cluster_info<-DT::renderDataTable({
    refresh$dynamic
    Cluster_id = unique(cyto@ClusterID)
    ID2CL=cyto@ID2CL
    inputLabel = c()
    for(i in 1:length(Cluster_id[,1])){
      inputLabel =c(inputLabel,as.character(ID2CL[which(ID2CL[,1]==Cluster_id[i,1]),2]))
    }
    size = c()
    for(i in 1:length(Cluster_id[,1])){
      size = c(size,length(cyto@label[which(cyto@label == as.character(Cluster_id[i,])),]))
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

  ## cell_type_confirm_plot ##
  output$cell_type_confirm_plot<-renderPlot({
    refresh$dynamic
    clusters=unlist(cyto@label)
    unique_cluster=unique(clusters)
    clusters[which(clusters!=unique_cluster[input$Cluster_info_rows_selected])]="Others"
    data=cyto@fcs.data[cyto@event.use,][,cyto@channel.use]

    Dimension1<-data[,input$confirm_marker_x]
    Dimension2<-data[,input$confirm_marker_y]

    if(is.null(input$confirm_marker_x)){
      return()
    }

    ggplot(data = as.data.frame(data),aes(x=Dimension1,y=Dimension2))+geom_point(aes(colour=clusters))+theme_bw()+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())+
      theme(panel.border  = element_blank())+
      theme(axis.line = element_line(colour = "black"))+
      xlab(input$confirm_marker_x)+
      ylab(input$confirm_marker_y)+
      scale_colour_manual(values = c("#3686ff","#d2d2d2"))
  })


  ## cell_type_confirm_plot_marker choose ##
  output$confirm_marker_choose<-renderUI({
    refresh$dynamic
    fluidRow(
      column(
        width = 6,
        selectInput(inputId = "confirm_marker_x",
                    label = "X",
                    selected = colnames(cyto@fcs.data[cyto@event.use,])[cyto@channel.use[1]],
                    choices = colnames(cyto@fcs.data[cyto@event.use,])[cyto@channel.use])
             ),
      column(
        width = 6,
        selectInput(inputId = "confirm_marker_y",
                    label = "Y",
                    selected = colnames(cyto@fcs.data[cyto@event.use,])[cyto@channel.use[2]],
                    choices = colnames(cyto@fcs.data[cyto@event.use,])[cyto@channel.use])
      )
    )
  })

  ##### confirm button #####
  observeEvent(input$confirm,{
    refresh$dynamic=refresh$dynamic+0.01
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

  #####=================== DMT plot page ===================#####
  output$DMT<-renderVisNetwork({
    refresh$dynamic
    if(is.null(cyto@dmt)){
      return()
    }

    #vector for nodes
    id=c(0)
    label=c("Root")
    title=c("Root")
    size=c(length(cyto@fcs.data[cyto@event.use,][,1]))
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

  output$cell_type_marker_line<-renderEChart({
    refresh$dynamic
    if(input$IDexchangeLabel=="ID"){
      ClusterID=cyto@ClusterID
    }
    else if(input$IDexchangeLabel=="Cell_Label"){
      ClusterID=cyto@label
    }

    cl_average_vec<-getAvg(fcsData=cyto@fcs.data[cyto@event.use,][,cyto@channel.use],ClusterID=ClusterID)
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
      setLegend(pos=3) %>%
      setTheme('macarons') %>%
      setToolbox(language = "en") %>%
      setGrid(width = 600,height = 580)
  })

 #####=================== MST plot page ===================#####
  output$MST_pie<-renderPlot({
    refresh$dynamic
    if(input$IDexchangeLabel=="ID"){
        ClusterID=cyto@ClusterID
      }
      else if(input$IDexchangeLabel=="Cell_Label"){
        ClusterID=cyto@label
      }

      MST=cyto@mst
        PlotPies(fsom = MST,cellTypes = as.matrix(ClusterID),view = "MST",
                 main = "Cluster pie MST Plot",
                 colorPalette = grDevices::colorRampPalette(c("#8DD3C7",
                                                              "#FFFFB3",
                                                              "#BEBADA",
                                                              "#FB8072",
                                                              "#80B1D3",
                                                              "#FDB462",
                                                              "#B3DE69",
                                                              "#FCCDE5",
                                                              "#D9D9D9",
                                                              "#BC80BD",
                                                              "#CCEBC5",
                                                              "#FFED6F")))
  })

  output$MST_marker<-renderPlot({
    refresh$dynamic
    MST=cyto@mst
    PlotMarker(MST,marker = input$MST_markers,view = "MST",
               main="Marker MST Plot",
               colorPalette = grDevices::colorRampPalette(c("#dbdcd7",
                                                            "#dddcd7",
                                                            "#e2c9cc",
                                                            "#e7627d",
                                                            "#b8235a",
                                                            "#801357",
                                                            "#3d1635",
                                                            "#1c1a27")))
  })


 ### MST pannel ###
  output$MST_panel<-renderUI({
    refresh$dynamic
    tagList(
      selectInput(inputId = "MST_markers",
                  label = "Markers",
                  choices = as.matrix(colnames(cyto@fcs.data[cyto@event.use,][,cyto@channel.use])),
                  selected = colnames(cyto@fcs.data[cyto@event.use,][,cyto@channel.use])[1])
    )
  })

 #####=================== Marker heatmap page ===================#####
  output$Marker_heatmap<-renderPlot({
    refresh$dynamic
    if(input$IDexchangeLabel=="ID"){
      ClusterID=cyto@ClusterID
    }
    else if(input$IDexchangeLabel=="Cell_Label"){
      ClusterID=cyto@label
    }
    cluster_av=getAvg(cyto@fcs.data[cyto@event.use,][,cyto@channel.use],ClusterID)
	if(input$heatmap_colorsets=="set1"){
		set=c("#ffffd9","#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58")
	}
	else{
		set=c("#1d56f8","#7379f4","#fefefe","#ef4c44","#f1345c")
	}
    pheatmap::pheatmap(cluster_av,
						color = colorRampPalette(set)(100),
						clustering_distance_rows=input$heatmap_clustering_distance_rows,
						clustering_distance_cols=input$heatmap_clustering_distance_cols,
						scale=input$heatmap_scale,
						fontsize=input$heatmap_fontsize
	)
  })


 #####=================== scatter plot page  ===================#####
  output$scatter_plot<-renderPlot({
    refresh$dynamic
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
    Expression = cyto@fcs.data[cyto@event.use,][input$Marker_exp_sca_choice]
    Dimension1=data[,1]
    Dimension2=data[,2]
    if(input$Cluster_Marker=="Markers"){
      set=c()
      if(input$scatter_colors=="set1"){
        set=c("#ffffd9","#edf8b1","#c7e9b4","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#253494","#081d58")
      }
      else if(input$scatter_colors=="set2"){
        set=c("#1d56f8","#7379f4","#fefefe","#ef4c44","#f1345c")
      }

      ggplot(data = data,aes(x=Dimension1,y=Dimension2))+geom_point(aes(colour=Expression))+theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())+
        theme(panel.border  = element_blank())+
        theme(axis.line = element_line(colour = "black"))+
        scale_color_gradientn(colours = set)
    }
    else if(input$Cluster_Marker=="Clusters"){
      if(is.null(input$Cluster_sca_choice)){
        return()
      }
      set=c()
      if(input$scatter_colors=="set1"){
        set=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")
      }
      else if(input$scatter_colors=="set2"){
        set=c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2")
      }
      if(input$Cluster_sca_choice=="ALL"){
        if(input$IDexchangeLabel=="ID"){
            clusters=as.factor(unlist(cyto@ClusterID))
        }
        else if(input$IDexchangeLabel=="Cell_Label"){
            clusters=as.factor(unlist(cyto@label))
        }

        ggplot(data = data,aes(x=Dimension1,y=Dimension2,colour=clusters))+geom_point()+theme_bw()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
          theme(panel.border = element_blank())+
          theme(axis.line = element_line(colour = "black"))+
          scale_colour_manual(values=colorRampPalette(set)(length(unique(clusters))))
      }
      else{
        if(input$IDexchangeLabel=="ID"){
            clusters=as.factor(unlist(cyto@ClusterID))
        }
        else if(input$IDexchangeLabel=="Cell_Label"){
            clusters=as.factor(unlist(cyto@label))
        }
        clusters=as.matrix(clusters)
        clusters[which(clusters!=input$Cluster_sca_choice)]="Others"
        clusters=as.factor(unlist(clusters))
        ggplot(data = data,aes(x=Dimension1,y=Dimension2,colour=clusters))+geom_point()+theme_bw()+
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank())+
          theme(panel.border  = element_blank())+
          theme(axis.line = element_line(colour = "black"))+
          scale_colour_manual(values = c("#3686ff","#d2d2d2"))
      }
    }
  })

 ###  markers or clusters event###
 observeEvent(input$Cluster_Marker,{
   refresh$dynamic
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
   refresh$dynamic
   if(input$IDexchangeLabel=="ID"){
     clusters=as.character(sort(unique(cyto@ClusterID)[,1]))
   }
   else if(input$IDexchangeLabel=="Cell_Label"){
     clusters=as.character(sort(unique(cyto@label)[,1]))
   }

  tagList(
    div(
      id="Cluster_sca",
      selectInput(inputId = "Cluster_sca_choice",
                  label = "you can select a cluster for displaying:",
                  choices = c("ALL",as.character(sort(unique(clusters)))),
                  selected = "ALL")
    ),
    div(
      id="Marker_exp_sca",style="display:none",
      selectInput(inputId = "Marker_exp_sca_choice",
                  label = "you can select a Marker for displaying:",
                  choices = c(sort(colnames(cyto@fcs.data[cyto@event.use,][,cyto@channel.use]))),
                  selected=cyto@fcs.data[cyto@event.use,][,cyto@channel.use][1])
      )
    )
  })

 #####=================== report page ===================#####

 # project part
 output$Project_report<-DT::renderDataTable({
   refresh$dynamic
   projectName <- cyto@projectname
   filename <- cyto@filename
   cell_num<-length(cyto@fcs.data[,1])
   cluster_num<-length(unique(cyto@label)[,1])
   transform_method<-cyto@transform_method
   channels<-paste(colnames(cyto@fcs.data[cyto@event.use,])[cyto@channel.use],collapse = ",\t")
   eventused<-length(cyto@event.use)
   clust_method<-cyto@clust_method
   preprocess<-cyto@preprocess
   AutoLabel<-cyto@autoLabel
   table<-data.frame(c("ProjectName","FileName","Cell num","num of clusters","Channels","events used","Transform method","Clust Method","Preprocess Method","AutoLabel"),
                     c(projectName,filename,cell_num,cluster_num,channels,eventused,transform_method,clust_method,preprocess,AutoLabel)
   )
   colnames(table)=c("Characters","Values")

   DT::datatable(
     table,
     class = "hover",
     style = "default",
     rownames=FALSE,
     escape = FALSE,
     selection="none",
     options = list(dom="lpt",paging = FALSE,ordering=FALSE,searching=FALSE,scrollY=400)
   )
 })
 # clusters part
 output$Cluster_report<-DT::renderDataTable({
   refresh$dynamic
   if(length(cyto@CL_label)==0){return()}
   refresh$dynamic
   ClusterID = unique(cyto@ClusterID)
   ID2CL=cyto@ID2CL
   ID = c()
   size = c()
   label = c()
   DMT_PhenoType=c()
   phenolist=list()
   for(i in 1:length(cyto@CL_label)){
     phenolist[cyto@CL_label[[i]]$ClusterID]=cyto@CL_label[[i]]$PhenoType
   }
   for(i in 1:length(ClusterID[,1])){
     ID = c(ID,ClusterID[i,1])
     size = c(size,length(which(cyto@ClusterID == ClusterID[i,1])))
     label =c(label,as.character(ID2CL[which(ID2CL[,1]==ClusterID[i,1]),2]))
     DMT_PhenoType=c(DMT_PhenoType,as.character(phenolist[as.character(ClusterID[i,1])]))
   }
   table=data.frame(ClusterID,size,label,DMT_PhenoType)

   data.frame(table)
 }, selection = "none",
   class = "hover",
   server = FALSE,
   escape = FALSE,
   rownames = FALSE,
   style = "default",
   options = list(
    dom = "lpt",
    paging = FALSE,
    scrollY=400
  )
 )

 # choose display report table
 observeEvent(input$choose_report,{
   if(input$choose_report=="Project"){
     showElement("Project_report_div")
     hideElement("Cluster_report_div")
   }
   else{
     showElement("Cluster_report_div")
     hideElement("Project_report_div")
   }
 })




 ## button for re-labeling the cell type ##
 observeEvent(input$back_labeling,{
   shinyalert(
     title = "GO back?",
     text = "Click 'Confirm' will go back to the labeling step and data produced from labeling will all be deleted from cytosee object",
     type = "warning",
     confirmButtonCol = "#3c8dbc",
     confirmButtonText = "Confirm",
     showCancelButton = TRUE,
     callbackR = function(x) { if(x != FALSE){
       cyto@label=list()
       Cluster_id=cyto@ClusterID
       cyto@ID2CL=cyto@ID2CL<<-data.frame(ID=Cluster_id,CD_Label = Cluster_id)
       showElement(id="change_marker_name")
       hideElement(id="celltype_tree_div")
       updateTabItems(session,"result_display",selected="ClusterLabel")
      }
     }
   )
 })


 ##### download Handler #######
 output$DownloadRdata <- downloadHandler(
   filename = function(){
     paste(cyto@projectname, ".RData", sep="")
   },
   content = function(file){
     refresh$dynamic
     save(cyto,file = file)
   }
 )

 #output$Img <- downloadHandler()
 }
   shinyApp(ui=ui, server = server,options = list("launch.browser"=TRUE))
}
