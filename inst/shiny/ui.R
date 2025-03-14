library(shinythemes)
library(shinyjs)


shinyUI(fluidPage(theme=shinytheme("cosmo"),
                  tags$head(tags$script(src="ExchangeData.js")),
                  useShinyjs(),  # Include shinyjs
                  extendShinyjs(script="ExchangeData.js", functions=c("send_results")),
                  tags$head(tags$style("body {background-image: url('/BackgroundTexture.jpg');background-size:cover;} 
  .shiny-myframe{padding:20px}
  .shiny-input-panel{background-color: transparent;border-color: transparent}
  .nav-tabs,.shiny-myframe{border-radius:5px;background-color:#E6E6EF;box-shadow: 2px 2px 4px #999999;padding:5px;}")),
                  
                  
                  h1(img(src="Logo.svg",style="width:150px"),"Variance-sensitive fuzzy clustering",style="text-shadow: 1px 0px #999900;font-weight:bold;color: #888888;"),
                p(paste0("Version ",packageVersion("vsclust")),br(), "See also: ",
                  a(href="https://github.com/computproteomics/vsclust","Source code and installation"),br(),"Publication: ",
                  a(href="https://doi.org/10.1093/bioinformatics/bty224","Veit Schwämmle, Ole N Jensen; Bioinformatics 2018,bty224"),br(),"Tutorial: ",a(href="https://pubmed.ncbi.nlm.nih.gov/33950508/", "Quantitative Methods in Proteomics"),style="font-size:8px;text-align:right"),
                  tabsetPanel(
                    tabPanel("File input", br(), 
                             p(
                               h2("File input"),
                               fileInput("in_file","Input file:",accept=c("txt/csv", "text/comma-separated-values,text/plain",".csv")),
                               textOutput("fileInText"),br(),
                               actionLink("reset", "Trigger server to reset file input"),br(),
                               actionLink("examplefile","load example"),
                               checkboxInput(inputId="is_header", label="Column names?", value=TRUE),
                               checkboxInput(inputId="protnames", label="Gene/protein identifiers in second column?", value=FALSE),
                               
                               #       sliderInput("fuzzifier",min=1.001,max=5,value=2,label="Fuzzifier value",step=0.001),
                               checkboxInput(inputId="isStat", label="Estimate variance levels from replicated quantifications? Otherwise: file contains mean values and one variance estimate per feature.",value=T)),hr(),
                             p(        h2("Experimental setup"),
                                       checkboxInput(inputId="isPaired", label="Paired tests",value=F),
                                       textOutput("RepsCond"),
                                       checkboxInput(inputId="qcol_order", label="Replicates are grouped",value = T),
                                       numericInput("NumReps",min=2,max=20,value=2,label="Number of replicates",step=1),
                                       numericInput("NumCond",min=2,max=20,value=3,label="Number of conditions",step=1)
                                       
                             ),hr(),value="fin"),
                    tabPanel("Statistics and variance",br(),htmlOutput("data_summ"),br(),
                             div(plotOutput("plot0",height=800),class="shiny-myframe"),br(),
                             downloadButton('downloadDataLimma', 'Download q-values and mean log-values '),
                             checkboxInput(inputId="limma_tsv", label="Download as tab-delimited file (default is a comma-separated file, csv)",value = F),
                             value = "stat"),
                    tabPanel("Estimation of cluster number",id="test",br(),inputPanel(
                      p("Maximum number of clusters for estimation:"),
                      sliderInput("maxclust",min=3,max=40,value=20,label=NULL, step=1),
                      actionButton("clButton1","Estimate parameters")),hr(),
                      div(plotOutput("plot1",height=600),class="shiny-myframe"), 
                      downloadButton('downloadParamEst', 'Download cluster figure'),br(),
                      value="pest"),
                    
                    tabPanel("Clustering results / variance-based", inputPanel(
                      p("Number of clusters"),
                      sliderInput("nclust1",min=3,max=40,value=5,label=NULL,step=1),
                      actionButton("clButton2","Run clustering")),hr(),
                      div(plotOutput("plot2"),class="shiny-plot-output"),
                      br(),
                      downloadButton('downloadData2', 'Download results '),
                      downloadButton('downloadCentroid2', 'Download centroids '),
                      downloadButton('downloadFigure', 'Download cluster figure'),
                      checkboxInput(inputId="clustvar_tsv", label="Download as tab-delimited file (default is a comma-separated file, csv)",value = F),
                      img(src="colormap.png",width="100px",height="50px"),
                      br(),
                      h3("Distribution of features over clusters"),
                      DT::dataTableOutput("clustinf1"),
                      value="clust1"),
                    tabPanel("Clustering results / standard method",inputPanel(
                      p("Number of clusters"),
                      sliderInput("nclust2",min=3,max=40,value=5,label=NULL,step=1),
                      actionButton("clButton3","Run clustering")),hr(),
                      div(plotOutput("plot3"),class="shiny-myframe"),
                      br(),
                      checkboxInput(inputId="cluststd_tsv", label="Download as tab-delimited file (default is a comma-separated file, csv)",value = F),
                      downloadButton('downloadData3', 'Download results'),
                      downloadButton('downloadCentroid3', 'Download centroids '),
                      downloadButton('downloadFigure2', 'Download cluster figure'),
                      img(src="colormap.png",width="100px",height="50px"),
                      br(),
                      h3("Distribution of features over clusters"),
                      DT::dataTableOutput("clustinf2"),
                      value="clust2"),
                    tabPanel("Enriched terms (DAVID)",inputPanel(
                      # radioButtons("enrich_method","Enrichment tool",choices=c("DAVID"="DAVID","KEGG (local copy)"="KEGG","GO molecular function (local copy)"="GOMF"),selected = "DAVID"),br(),
                      selectInput("infosource","Information resource (DAVID)",choices=list("GO terms"=c("GO molecular function"="Function",
                                                                                                        "GO biological process"="Process",
                                                                                                        "GO cellular component"="Component"),
                                                                                           UniProt = c("UniProt keywords" = "Keyword"),
                                                                                           Pathways=c("KEGG"="KEGG",
                                                                                                      # "PANTHER"="PANTHER_PATHWAY",
                                                                                                      "REACTOME"="RCTM",
                                                                                                      "Wikipathways"="Wikipathways"
                                                                                           ),
                                                                                           "Protein domains"=c("Pfam"="Pfam","SMART"="SMART",
                                                                                                                   "InterPro"="InterPro","BioGRID"="BIOGRID_INTERACTION"),
                                                                                           "TISSUES" = c("TISSUES"),
                                                                                           "DISEASES" = c("DISEASES"),
                                                                                           "COMPARTMENTS" = c("COMPARTMENTS")
                                                                                           
                      ), multiple = T,selected="Process"),

                      # selectInput("organism","Organism (not needed for DAVID)",c("anopheles","arabidopsis","bovine","canine", "chicken", "chimp", "coelicolor", "ecolik12","ecsakai", "fly", "gondii","human", "malaria", "mouse", "pig", "rat","rhesus", "worm", "xenopus", "yeast","zebrafish"),selected="human"),
                      actionButton("goButton","Run enrichment")),hr(),
                      h3("Variance-based clustering"),
                      div(plotOutput("plot4"),class="shiny-myframe"),
                      p("Plot shows only top 20 terms"),
                      downloadButton('downloadGOData1', 'Download full results\n(variance-based clustering)'),br(),
                      h3("Standard clustering"),
                      div(plotOutput("plot5"),class="shiny-myframe"),
                      p("Plot shows only top 20 terms"),          
                      downloadButton('downloadGOData2', 'Download full results\n(standard clustering)'),
                      value="gos"),
                    tabPanel("Help",h4("Introduction"),
                             htmlOutput("intro"),
                             downloadLink("DownloadExample","Download example file"),
                             h4("File input"),htmlOutput("finput"),
                             h4("Statistical analysis"),htmlOutput("stat"),
                             h4("Parameter estimation"),htmlOutput("pest"),
                             h4("Clustering"),htmlOutput("fclust"),
                             h4("GO terms"),htmlOutput("goterms"),
                             h4("Further information"),htmlOutput("reading"), value="help"),id="tabset")
)

)
