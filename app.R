# app.R
library(shiny)
library(lme4)
library(emmeans)
library(DT) # Replaces the deprecated Shiny tables

# --- USER INTERFACE ---
ui <- bootstrapPage(
  
  # THIS INJECTS TAILWIND CSS DIRECTLY!
  tags$head(tags$script(src = "https://cdn.tailwindcss.com")),
  
  div(class = "min-h-screen bg-gray-50 p-4 md:p-8 font-sans",
      
      div(class = "max-w-6xl mx-auto space-y-6",
          
          # Header Card
          div(class = "bg-white rounded-xl shadow-sm p-6 border border-gray-100",
              h1(class = "text-2xl font-bold text-blue-700 mb-2", "Multi-Environment Trial Analyzer (Two-Stage Model)"),
              p(class = "text-gray-500 text-sm", "Note: Processing complex LMMs in the browser may take a few moments depending on dataset size.")
          ),
          
          # Upload Controls
          div(class = "bg-white rounded-xl shadow-sm p-6 border border-gray-100 flex flex-col md:flex-row gap-4 items-end",
              div(class = "flex-1 w-full",
                  fileInput("file", "Upload Trial Data (CSV)", accept = ".csv", width = "100%")
              ),
              div(class = "w-full md:w-auto mb-4",
                  actionButton("run", "Run Two-Stage Analysis", class = "w-full bg-blue-600 hover:bg-blue-700 text-white font-semibold py-2 px-6 rounded-lg shadow transition-colors")
              )
          ),
          
          # Main panel for results
          uiOutput("error_ui"),
          uiOutput("results_ui")
      )
  )
)

# --- SERVER LOGIC ---
server <- function(input, output, session) {
  
  analysis_results <- reactiveVal(NULL)
  error_msg <- reactiveVal(NULL)
  
  # Trigger analysis on button click
  observeEvent(input$run, {
    req(input$file)
    error_msg(NULL)
    analysis_results(NULL)
    
    tryCatch({
      data <- read.csv(input$file$datapath)
      
      # Factorize columns
      cols_to_factor <- c("Name", "Location", "Season", "Block", "Rep")
      has_trial <- "Trial.Code" %in% names(data)
      if(has_trial) cols_to_factor <- c(cols_to_factor, "Trial.Code")
      
      cols_to_factor <- intersect(cols_to_factor, names(data))
      data[cols_to_factor] <- lapply(data[cols_to_factor], as.factor)
      data$YLD_MKT.P <- as.numeric(data$YLD_MKT.P)
      
      if(!has_trial) {
        data$Trial.Code <- as.factor(paste(data$Location, data$Season, sep="_"))
      }
      
      # STAGE 1
      split_data <- split(data, data$Trial.Code)
      stage1_results <- list()
      
      for(trial_name in names(split_data)) {
        d <- split_data[[trial_name]]
        if(length(unique(d$Name)) < 2) next 
        
        tryCatch({
          if("Block" %in% names(d) && length(unique(d$Block)) > 1) {
            m1 <- lmer(YLD_MKT.P ~ Name + (1|Rep) + (1|Rep:Block), data = d)
          } else {
            m1 <- lmer(YLD_MKT.P ~ Name + (1|Rep), data = d)
          }
          
          em <- as.data.frame(emmeans(m1, "Name"))
          em$Weight <- 1 / (em$SE^2)
          
          stage1_results[[trial_name]] <- data.frame(
            Trial_Code = trial_name,
            Genotype = em$Name,
            BLUE = round(em$emmean, 3),
            SE = round(em$SE, 3),
            Weight = em$Weight
          )
        }, error = function(e) {
          # Skip singular fits silently for cleaner UI
        })
      }
      
      stage1_df <- do.call(rbind, stage1_results)
      rownames(stage1_df) <- NULL
      
      # STAGE 2
      grand_mean <- mean(stage1_df$BLUE, na.rm = TRUE)
      m2 <- lmer(BLUE ~ (1|Genotype) + (1|Trial_Code), weights = Weight, data = stage1_df)
      
      ranef_m2 <- ranef(m2, condVar = TRUE)
      blup_values <- ranef_m2$Genotype
      
      post_vars <- attr(blup_values, "postVar")
      seps <- sqrt(post_vars[1, 1, ])
      
      unadj_means <- aggregate(YLD_MKT.P ~ Name, data = data, FUN = function(x) round(mean(x, na.rm=TRUE), 3))
      names(unadj_means) <- c("Genotype", "Unadj_Mean")
      
      blup_df <- data.frame(
        Genotype = rownames(blup_values),
        BLUP_Adj = round(blup_values$`(Intercept)`, 3),
        SEP = round(seps, 3), 
        Predicted_Yield = round(grand_mean + blup_values$`(Intercept)`, 3)
      )
      
      blup_df <- merge(blup_df, unadj_means, by = "Genotype")
      blup_df$Raw_Deviation <- round(blup_df$Unadj_Mean - grand_mean, 3)
      blup_df <- blup_df[order(-blup_df$Predicted_Yield), ]
      rownames(blup_df) <- NULL
      
      var_comp <- as.data.frame(VarCorr(m2))
      var_g <- round(var_comp$vcov[var_comp$grp == "Genotype"], 4)
      var_trial <- round(var_comp$vcov[var_comp$grp == "Trial_Code"], 4)
      var_err <- round(var_comp$vcov[var_comp$grp == "Residual"], 4)
      
      cullis_h2 <- if (var_g > 0) round(1 - (mean(seps^2) / var_g), 3) else 0
      
      # Save to reactive
      analysis_results(list(
        stats = list(GM = round(grand_mean, 2), H2 = cullis_h2, Vg = var_g, Vt = var_trial, Ve = var_err),
        s1 = stage1_df,
        s2 = blup_df
      ))
      
    }, error = function(e) {
      error_msg(e$message)
    })
  })
  
  # --- UI RENDERING ---
  output$error_ui <- renderUI({
    req(error_msg())
    div(class = "bg-red-50 border border-red-200 text-red-700 px-4 py-3 rounded-lg mt-4", error_msg())
  })
  
  output$results_ui <- renderUI({
    req(analysis_results())
    res <- analysis_results()
    
    tagList(
      # Top level stats
      div(class = "bg-white rounded-xl shadow-sm border border-blue-200 mb-8 overflow-hidden",
          div(class = "bg-blue-600 text-white px-6 py-3 font-semibold", "1. Overall Multi-Environment Variance (Stage 2)"),
          div(class = "p-6",
              div(class = "grid grid-cols-1 md:grid-cols-3 gap-4 mb-6",
                  div(class = "bg-blue-50 rounded-lg p-4 text-center border border-blue-100",
                      div(class = "text-xs text-blue-600 font-bold uppercase tracking-wider", "Genotype Variance"),
                      div(class = "text-3xl font-bold text-blue-900 mt-1", res$stats$Vg)
                  ),
                  div(class = "bg-gray-50 rounded-lg p-4 text-center border border-gray-200",
                      div(class = "text-xs text-gray-500 font-bold uppercase tracking-wider", "Trial Variance"),
                      div(class = "text-3xl font-bold text-gray-800 mt-1", res$stats$Vt)
                  ),
                  div(class = "bg-indigo-50 rounded-lg p-4 text-center border border-indigo-100",
                      div(class = "text-xs text-indigo-600 font-bold uppercase tracking-wider", "Residual Variance"),
                      div(class = "text-3xl font-bold text-indigo-900 mt-1", res$stats$Ve)
                  )
              ),
              hr(class = "border-gray-200 mb-4"),
              h4(class = "text-lg text-gray-800", HTML(paste0("Overall Cullis Heritability (H&sup2;): <strong class='text-blue-600'>", res$stats$H2, "</strong>")))
          )
      ),
      
      # Stage 1 Table
      div(class = "bg-white rounded-xl shadow-sm border border-gray-200 mb-8 overflow-hidden",
          div(class = "bg-gray-50 px-6 py-3 border-b border-gray-200 flex justify-between items-center",
              span(class = "font-semibold text-gray-800", "2. Stage 1: Trial BLUEs & Reliability Weights"),
              downloadButton("dl_s1", "Download BLUEs", class = "bg-white hover:bg-gray-100 text-blue-600 font-medium py-1 px-3 border border-blue-200 rounded text-sm transition")
          ),
          div(class = "p-6",
              p(class = "text-sm text-gray-500 mb-4", "Notice the Weight column. Trials with high Standard Errors (SE) receive tiny weights."),
              DTOutput("tbl_s1")
          )
      ),
      
      # Stage 2 Table
      div(class = "bg-white rounded-xl shadow-sm border border-gray-200 mb-8 overflow-hidden",
          div(class = "bg-gray-50 px-6 py-3 border-b border-gray-200 flex justify-between items-center",
              span(class = "font-semibold text-gray-800", "3. Final Multi-Environment BLUPs (Stage 2)"),
              downloadButton("dl_s2", "Download BLUPs", class = "bg-white hover:bg-gray-100 text-green-600 font-medium py-1 px-3 border border-green-200 rounded text-sm transition")
          ),
          div(class = "p-6",
              p(class = "text-sm text-gray-500 mb-4", HTML(paste0("Grand Mean across all weighted trials is <strong class='text-gray-900'>", res$stats$GM, "</strong>."))),
              DTOutput("tbl_s2")
          )
      )
    )
  })
  
  # Table Renderers
  output$tbl_s1 <- renderDT({
    req(analysis_results())
    analysis_results()$s1
  }, options = list(
    pageLength = 50, 
    scrollX = TRUE,
    scrollY = "400px",      
    scrollCollapse = TRUE   
  ))
  
  output$tbl_s2 <- renderDT({
    req(analysis_results())
    analysis_results()$s2
  }, options = list(
    pageLength = 50, 
    scrollX = TRUE          
  ))
  
  # Download Handlers
  output$dl_s1 <- downloadHandler(
    filename = function() { "stage1_trial_blues_weights.csv" },
    content = function(file) { write.csv(analysis_results()$s1, file, row.names = FALSE) }
  )
  
  output$dl_s2 <- downloadHandler(
    filename = function() { "overall_twostage_blups.csv" },
    content = function(file) { write.csv(analysis_results()$s2, file, row.names = FALSE) }
  )
}

shinyApp(ui, server)