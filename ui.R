# Load required libraries
library(shiny)
library(shinydashboard)
library(DT)

# =============================================================================
# UI Definition (保持原有UI不变)
# =============================================================================

ui <- dashboardPage(
  dashboardHeader(title = "scGTM: Single-cell Gene Trend Model with Cuckoo Search"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Gene Analysis", tabName = "analysis", icon = icon("dna")),
      menuItem("Biological Insights", tabName = "insights", icon = icon("microscope")),
      menuItem("Batch Comparison", tabName = "batch", icon = icon("layer-group")),
      menuItem("Help & Guide", tabName = "help", icon = icon("question"))
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .main-header .navbar { background-color: #2c3e50 !important; }
        .content-wrapper { background-color: #ecf0f1; }
        .box.box-primary { border-top-color: #3498db; }
        .box.box-success { border-top-color: #27ae60; }
        .box.box-warning { border-top-color: #f39c12; }
        .box.box-info { border-top-color: #3498db; }
        .tip-box {
          background: #e8f8f5;
          border-left: 4px solid #27ae60;
          padding: 10px;
          margin: 10px 0;
        }
        .warning-box {
          background: #fef9e7;
          border-left: 4px solid #f39c12;
          padding: 10px;
          margin: 10px 0;
        }
      "))
    ),
    
    tabItems(
      # Data Upload Tab
      tabItem(tabName = "upload",
              fluidRow(
                box(title = "Data Loading", status = "primary", solidHeader = TRUE, width = 12,
                    column(6,
                           fileInput("data_file", "Upload CSV File", accept = c(".csv")),
                           checkboxInput("header", "Header", TRUE),
                           checkboxInput("use_example", "Use Example Data", FALSE),
                           br(),
                           actionButton("load_data", "Load Data", class = "btn-primary"),
                           br(), br(),
                           div(class = "tip-box",
                               h5(icon("info"), " Data Format Requirements:"),
                               tags$ul(
                                 tags$li("Column 1: Cell identifiers"),
                                 tags$li("Column 2: Pseudotime values (0-1)"),
                                 tags$li("Columns 3+: Gene expression counts")
                               )
                           )
                    ),
                    column(6,
                           h4("Data Summary"),
                           verbatimTextOutput("data_summary")
                    )
                )
              ),
              
              fluidRow(
                box(title = "Data Preview", status = "info", width = 12,
                    DT::dataTableOutput("data_preview")
                )
              )
      ),
      
      # Gene Analysis Tab
      tabItem(tabName = "analysis",
              fluidRow(
                box(title = "Gene Selection & Model Configuration", status = "primary", 
                    solidHeader = TRUE, width = 4,
                    selectInput("gene_select", "Select Gene for Analysis:",
                                choices = NULL),
                    
                    h4("Model Selection"),
                    checkboxGroupInput("selected_models", "Select models to compare:",
                                       choices = list(
                                         "Hill-shaped (peak)" = "hill",
                                         "Valley-shaped (dip)" = "valley", 
                                         "Monotonic Increasing" = "increasing",
                                         "Monotonic Decreasing" = "decreasing"
                                       ),
                                       selected = c("hill", "valley", "increasing", "decreasing")),
                    
                    h4("Distribution Selection"),
                    selectInput("distribution", "Select count distribution:",
                                choices = list(
                                  "Poisson" = "poisson",
                                  "Negative Binomial" = "nb",
                                  "Zero-Inflated Poisson" = "zip",
                                  "Zero-Inflated Negative Binomial" = "zinb"
                                ),
                                selected = "poisson"),
                    
                    div(class = "tip-box",
                        h5(icon("lightbulb"), " scGTM with Cuckoo Search:"),
                        tags$ul(
                          tags$li("Cuckoo Search optimization algorithm"),
                          tags$li("Multiple count distributions"),
                          tags$li("Automatic model selection via AIC"),
                          tags$li("Automatic asymptotic variance calculation")
                        )
                    ),
                    
                    h4("Cuckoo Search Parameters"),
                    numericInput("n_nests", "Number of Nests:", value = 25, min = 15, max = 50),
                    numericInput("max_iter", "Max Iterations:", value = 1000, min = 500, max = 2000),
                    numericInput("pa", "Discovery Rate (pa):", value = 0.25, min = 0.1, max = 0.5, step = 0.05),
                    
                    br(),
                    actionButton("fit_model", "Fit scGTM Model", 
                                 class = "btn-success btn-lg", width = "100%"),
                    br(), br(),
                    verbatimTextOutput("fitting_status")
                ),
                
                box(title = "Gene Expression Trend", status = "info", width = 8,
                    plotOutput("trend_plot", height = "400px"),
                    br(),
                    fluidRow(
                      column(6,
                             h5("Detected Pattern:"),
                             textOutput("detected_pattern")
                      ),
                      column(6,
                             h5("Distribution & Quality:"),
                             textOutput("distribution_info")
                      )
                    )
                )
              ),
              
              fluidRow(
                box(title = "Model Comparison", status = "info", width = 12,
                    h4("All Models Comparison:"),
                    tableOutput("model_comparison_table"),
                    br(),
                    h5("Model Selection Details:"),
                    verbatimTextOutput("model_selection_details")
                )
              ),
              
              fluidRow(
                box(title = "Model Parameters & Asymptotic Variance", status = "success", width = 6,
                    h4("Parameter Estimates:"),
                    tableOutput("parameter_table"),
                    br(),
                    h4("Asymptotic Standard Errors:"),
                    tableOutput("asymptotic_se_table")
                ),
                
                box(title = "Model Quality Assessment", status = "warning", width = 6,
                    verbatimTextOutput("model_quality"),
                    br(),
                    downloadButton("download_results", "Download Analysis", class = "btn-info")
                )
              ),
              
              fluidRow(
                box(title = "Model Diagnostics", status = "info", width = 12,
                    plotOutput("diagnostic_plots", height = "300px")
                )
              )
      ),
      
      # Biological Insights Tab
      tabItem(tabName = "insights",
              fluidRow(
                box(title = "Biological Pattern Classification", status = "primary", 
                    solidHeader = TRUE, width = 6,
                    verbatimTextOutput("pattern_classification"),
                    br(),
                    div(class = "tip-box",
                        h5(icon("info"), " Pattern Types from scGTM Paper:"),
                        tags$ul(
                          tags$li(strong("Hill-shaped:"), "Peak expression pattern"),
                          tags$li(strong("Valley-shaped:"), "Dip expression pattern"),
                          tags$li(strong("Monotonic Increasing:"), "Continuous upregulation"),
                          tags$li(strong("Monotonic Decreasing:"), "Continuous downregulation"),
                          tags$li(strong("Early-acting:"), "t₀ < 0.3"),
                          tags$li(strong("Late-acting:"), "t₀ > 0.7")
                        )
                    )
                ),
                
                box(title = "Expression Statistics", status = "info", width = 6,
                    tableOutput("expression_stats"),
                    br(),
                    plotOutput("pseudotime_dist", height = "200px")
                )
              ),
              
              fluidRow(
                box(title = "Mathematical Models", status = "success", width = 12,
                    div(style = "background: #f8f9fa; padding: 15px; border-radius: 5px;",
                        h5("Four scGTM Models:"),
                        div(style = "font-family: 'Courier New', monospace; font-size: 12px;",
                            strong("Hill:"), " log(E[Y]+1) = μ×exp(-k₁×(t-t₀)²) if t≤t₀, μ×exp(-k₂×(t-t₀)²) if t>t₀", br(),
                            strong("Valley:"), " log(E[Y]+1) = b-μ×exp(-k₂×(t-t₀)²) if t≤t₀, b-μ×exp(-k₁×(t-t₀)²) if t>t₀", br(),
                            strong("Increasing:"), " log(E[Y]+1) = μ×(1+tanh(k×(t-t₀)))/2", br(),
                            strong("Decreasing:"), " log(E[Y]+1) = μ×(1-tanh(k×(t-t₀)))/2"
                        )
                    ),
                    plotOutput("model_illustration", height = "300px")
                )
              )
      ),
      
      # Batch Comparison Tab
      tabItem(tabName = "batch",
              fluidRow(
                box(title = "Multi-Gene Comparison", status = "primary", 
                    solidHeader = TRUE, width = 12,
                    column(3,
                           h4("Gene Selection"),
                           checkboxGroupInput("batch_genes", "Select genes:",
                                              choices = NULL),
                           br(),
                           actionButton("select_all_genes", "Select All", class = "btn-info btn-sm"),
                           actionButton("clear_all_genes", "Clear All", class = "btn-warning btn-sm")
                    ),
                    
                    column(3,
                           h4("Algorithm Parameters"),
                           numericInput("batch_nests", "Nests per gene:", value = 20, min = 15, max = 35),
                           numericInput("batch_iterations", "Iterations per gene:", value = 800, min = 500, max = 1500),
                           selectInput("batch_distribution", "Distribution:",
                                       choices = list("Poisson" = "poisson", "NB" = "nb", "ZIP" = "zip", "ZINB" = "zinb"),
                                       selected = "poisson")
                    ),
                    
                    column(3,
                           h4("Model Selection"),
                           checkboxGroupInput("batch_models", "Models to compare:",
                                              choices = list(
                                                "Hill" = "hill",
                                                "Valley" = "valley",
                                                "Increasing" = "increasing",
                                                "Decreasing" = "decreasing"
                                              ),
                                              selected = c("hill", "valley"))
                    ),
                    
                    column(3,
                           h4("Control"),
                           br(),
                           actionButton("run_batch", "Run Analysis", class = "btn-success btn-lg"),
                           br(), br(),
                           downloadButton("download_batch", "Download Results", class = "btn-info"),
                           br(), br(),
                           verbatimTextOutput("batch_status")
                    )
                )
              ),
              
              fluidRow(
                box(title = "Comparative Results", status = "info", width = 12,
                    DT::dataTableOutput("batch_results_table"),
                    br(),
                    plotOutput("pattern_summary_plot", height = "250px")
                )
              ),
              
              fluidRow(
                box(title = "Gene Clustering by Expression Dynamics", status = "success", width = 12,
                    plotOutput("gene_clustering_plot", height = "400px")
                )
              )
      ),
      
      # Help Tab
      tabItem(tabName = "help",
              fluidRow(
                box(title = "About scGTM with Cuckoo Search", status = "primary", width = 6,
                    h3("What is scGTM?"),
                    p("The Single-cell Generalized Trend Model (scGTM) captures interpretable gene expression trends along cell pseudotime using four mathematical models and Cuckoo Search optimization."),
                    
                    h4("Key Features:"),
                    tags$ol(
                      tags$li("Four model types: Hill, Valley, Increasing, Decreasing"),
                      tags$li("Multiple count distributions: Poisson, NB, ZIP, ZINB"),
                      tags$li("Cuckoo Search optimization algorithm"),
                      tags$li("Automatic model selection via AIC"),
                      tags$li("Automatic asymptotic variance calculation using numDeriv")
                    ),
                    
                    h4("Cuckoo Search Algorithm:"),
                    p("Bio-inspired optimization algorithm that mimics cuckoo birds' brood parasitism behavior. Superior to PSO for complex optimization landscapes with automatic variance calculation.")
                ),
                
                box(title = "Usage Guide", status = "success", width = 6,
                    h4("Workflow:"),
                    tags$ol(
                      tags$li("Upload data (cells × genes + pseudotime)"),
                      tags$li("Select gene and models to compare"),
                      tags$li("Choose appropriate count distribution"),
                      tags$li("Configure Cuckoo Search parameters"),
                      tags$li("Fit scGTM model and interpret results")
                    ),
                    
                    h4("Distribution Selection Tips:"),
                    tags$ul(
                      tags$li(strong("Poisson:"), "Simple count data"),
                      tags$li(strong("Negative Binomial:"), "Overdispersed counts"),
                      tags$li(strong("ZIP:"), "Excess zeros in Poisson"),
                      tags$li(strong("ZINB:"), "Excess zeros + overdispersion")
                    )
                )
              ),
              
              fluidRow(
                box(title = "Parameter Interpretation", status = "info", width = 6,
                    h4("Model Parameters:"),
                    tags$ul(
                      tags$li(strong("μ:"), "Expression magnitude"),
                      tags$li(strong("k₁, k₂:"), "Steepness parameters"),
                      tags$li(strong("t₀:"), "Peak/valley/inflection time"),
                      tags$li(strong("k (monotonic):"), "Rate parameter")
                    ),
                    
                    h4("Quality Metrics:"),
                    tags$ul(
                      tags$li("R² > 0.8: Excellent fit"),
                      tags$li("R² > 0.6: Good fit"),
                      tags$li("AIC: Lower is better"),
                      tags$li("Asymptotic SE: Parameter uncertainty (auto-calculated)")
                    )
                ),
                
                box(title = "References", status = "primary", width = 6,
                    h4("Primary Reference:"),
                    div(style = "background: #f8f9fa; padding: 15px; border-radius: 5px;",
                        p(strong("Cui, E. H., Song, D., Wong, W. K., & Li, J. J. (2022)."), 
                          "Single-cell generalized trend model (scGTM): A flexible and interpretable model of gene expression trend along cell pseudotime. ", 
                          em("Bioinformatics"), ", 38(16), 3927-3934.")
                    ),
                    
                    h4("Implementation:"),
                    tags$ul(
                      tags$li("Cuckoo Search optimization"),
                      tags$li("Multiple count distributions"),
                      tags$li("numDeriv package for automatic Hessian calculation"),
                      tags$li("Automatic model selection"),
                      tags$li("R Shiny interface")
                    )
                )
              )
      )
    )
  )
)
