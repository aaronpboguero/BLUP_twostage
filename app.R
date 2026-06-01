# app.R
library(shiny)
library(lme4)
library(emmeans)
library(reactable)

# --- USER INTERFACE ---
ui <- bootstrapPage(
  
  tags$head(tags$script(src = "https://cdn.tailwindcss.com")),
  
  div(class = "min-h-screen bg-gray-50 p-4 md:p-8 font-sans",
      
      div(class = "max-w-6xl mx-auto space-y-6",
          
          # Header Card
          div(class = "bg-white rounded-xl shadow-sm p-6 border border-gray-100",
              h1(class = "text-2xl font-bold text-blue-700 mb-2", "Multi-Environment Trial Analyzer"),
              p(class = "text-gray-500 text-sm", "Upload your data. Include 'Female' and 'Male' columns to automatically extract GCA and SCA.")
          ),
          
          # Upload Controls
          div(class = "bg-white rounded-xl shadow-sm p-6 border border-gray-100 flex flex-col md:flex-row gap-4 items-end",
              div(class = "flex-1 w-full",
                  fileInput("file", "Upload Trial Data (CSV)", accept = ".csv", width = "100%")
              ),
              div(class = "w-full md:w-auto mb-4",
                  actionButton("run", "Run Analysis", class = "w-full bg-blue-600 hover:bg-blue-700 text-white font-semibold py-2 px-6 rounded-lg shadow transition-colors")
              )
          ),
          
          # Error display
          uiOutput("error_ui"),
          
          # DASHBOARD STRUCTURE
          conditionalPanel(
            condition = "output.analysis_ready",
            
            # 1. Stage 1 Table
            div(class = "bg-white rounded-xl shadow-sm border border-gray-200 mb-8 overflow-hidden",
                div(class = "bg-gray-50 px-6 py-3 border-b border-gray-200",
                    span(class = "font-semibold text-gray-800", "1. Stage 1: Trial BLUEs & Weights")
                ),
                div(class = "p-6", reactableOutput("tbl_s1"))
            ),
            
            # 2. Stage 2 Table (Overall BLUPs)
            div(class = "bg-white rounded-xl shadow-sm border border-gray-200 mb-8 overflow-hidden",
                div(class = "bg-gray-50 px-6 py-3 border-b border-gray-200",
                    span(class = "font-semibold text-gray-800", "2. Final Multi-Environment Hybrid BLUPs")
                ),
                div(class = "p-6", reactableOutput("tbl_s2"))
            ),
            
            # 3. Dynamic Parental/Combining Ability Tables
            uiOutput("ca_tables_ui")
          )
      )
  )
)

# --- SERVER LOGIC ---
server <- function(input, output, session) {
  
  analysis_results <- reactiveVal(NULL)
  error_msg <- reactiveVal(NULL)
  
  output$analysis_ready <- reactive({ !is.null(analysis_results()) })
  outputOptions(output, "analysis_ready", suspendWhenHidden = FALSE)
  
  observeEvent(input$run, {
    req(input$file)
    error_msg(NULL)
    analysis_results(NULL)
    
    tryCatch({
      data <- read.csv(input$file$datapath)
      
      # Determine if parents are present
      has_parents <- all(c("Female", "Male") %in% names(data))
      
      # Factorize columns
      cols_to_factor <- c("Name", "Location", "Season", "Block", "Rep", "Female", "Male")
      has_trial <- "Trial.Code" %in% names(data)
      if(has_trial) cols_to_factor <- c(cols_to_factor, "Trial.Code")
      
      cols_to_factor <- intersect(cols_to_factor, names(data))
      data[cols_to_factor] <- lapply(data[cols_to_factor], as.factor)
      data$YLD_MKT.P <- as.numeric(data$YLD_MKT.P)
      
      if(!has_trial) {
        data$Trial.Code <- as.factor(paste(data$Location, data$Season, sep="_"))
      }
      
      # --- STAGE 1: Calculate BLUEs per Trial ---
      split_data <- split(data, data$Trial.Code)
      stage1_results <- list()
      
      for(trial_name in names(split_data)) {
        d <- split_data[[trial_name]]
        if(length(unique(d$Name)) < 2) next 
        
        # Simple LMM for Stage 1
        m1 <- lmer(YLD_MKT.P ~ Name + (1|Rep), data = d)
        em <- as.data.frame(emmeans(m1, "Name"))
        raw_weight <- ifelse(is.na(em$SE) | em$SE == 0, 0, 1 / (em$SE^2)) 
        
        stage1_results[[trial_name]] <- data.frame(
          Trial_Code = trial_name,
          Genotype = em$Name,
          BLUE = round(em$emmean, 3),
          Weight = round(raw_weight, 4) 
        )
      }
      
      stage1_df <- do.call(rbind, stage1_results)
      rownames(stage1_df) <- NULL
      
      # Map parents back to Stage 1 data if they exist
      if(has_parents) {
        parent_map <- unique(data[, c("Name", "Female", "Male")])
        stage1_df <- merge(stage1_df, parent_map, by.x = "Genotype", by.y = "Name", all.x = TRUE)
      }
      
      # --- STAGE 2: Multi-Environment BLUPs ---
      grand_mean <- mean(stage1_df$BLUE, na.rm = TRUE)
      
      if(has_parents) {
        # NCII Combining Ability Model
        m2 <- lmer(BLUE ~ (1|Female) + (1|Male) + (1|Genotype) + (1|Trial_Code), weights = Weight, data = stage1_df)
        
        # Extract GCAs
        f_gca <- data.frame(Female = rownames(ranef(m2)$Female), GCA_Female = round(ranef(m2)$Female$`(Intercept)`, 3))
        f_gca <- f_gca[order(-f_gca$GCA_Female), , drop = FALSE]
        
        m_gca <- data.frame(Male = rownames(ranef(m2)$Male), GCA_Male = round(ranef(m2)$Male$`(Intercept)`, 3))
        m_gca <- m_gca[order(-m_gca$GCA_Male), , drop = FALSE]
        
        # Extract SCA
        sca_df <- data.frame(Hybrid = rownames(ranef(m2)$Genotype), SCA = round(ranef(m2)$Genotype$`(Intercept)`, 3))
        sca_df <- sca_df[order(-sca_df$SCA), , drop = FALSE]
        
      } else {
        # Standard Hybrid Model
        m2 <- lmer(BLUE ~ (1|Genotype) + (1|Trial_Code), weights = Weight, data = stage1_df)
        f_gca <- NULL
        m_gca <- NULL
        sca_df <- NULL
      }
      
      # Base Hybrid BLUP extraction
      blup_values <- ranef(m2)$Genotype
      blup_df <- data.frame(
        Genotype = rownames(blup_values),
        Predicted_Yield = round(grand_mean + blup_values$`(Intercept)`, 3)
      )
      blup_df <- blup_df[order(-blup_df$Predicted_Yield), ]
      rownames(blup_df) <- NULL
      
      # Save to reactive
      analysis_results(list(
        has_parents = has_parents,
        s1 = stage1_df,
        s2 = blup_df,
        f_gca = f_gca,
        m_gca = m_gca,
        sca = sca_df
      ))
      
    }, error = function(e) {
      error_msg(paste("Model Error:", e$message, "- Note: With only 1 loc/2 reps, lme4 may trigger a singular fit if variance is zero."))
    })
  })
  
  # --- UI RENDERING ---
  output$error_ui <- renderUI({
    req(error_msg())
    div(class = "bg-red-50 border border-red-200 text-red-700 px-4 py-3 rounded-lg mt-4", error_msg())
  })
  
  # Render Stage 1 Table
  output$tbl_s1 <- renderReactable({
    req(analysis_results())
    reactable(analysis_results()$s1, pagination = TRUE, defaultPageSize = 5, compact = TRUE)
  })
  
  # Render Stage 2 Table
  output$tbl_s2 <- renderReactable({
    req(analysis_results())
    reactable(analysis_results()$s2, pagination = TRUE, defaultPageSize = 5, compact = TRUE)
  })
  
  # Render Dynamic CA Tables
  output$ca_tables_ui <- renderUI({
    req(analysis_results())
    if(!analysis_results()$has_parents) return(NULL)
    
    div(class = "space-y-8",
        # Female GCA
        div(class = "bg-white rounded-xl shadow-sm border border-pink-200 overflow-hidden",
            div(class = "bg-pink-50 px-6 py-3 border-b border-pink-200",
                span(class = "font-semibold text-pink-800", "3A. Female Parent GCA")
            ),
            div(class = "p-6", reactableOutput("tbl_f_gca"))
        ),
        # Male GCA
        div(class = "bg-white rounded-xl shadow-sm border border-blue-200 overflow-hidden",
            div(class = "bg-blue-50 px-6 py-3 border-b border-blue-200",
                span(class = "font-semibold text-blue-800", "3B. Male Parent GCA")
            ),
            div(class = "p-6", reactableOutput("tbl_m_gca"))
        ),
        # SCA
        div(class = "bg-white rounded-xl shadow-sm border border-purple-200 overflow-hidden",
            div(class = "bg-purple-50 px-6 py-3 border-b border-purple-200",
                span(class = "font-semibold text-purple-800", "3C. Specific Combining Ability (SCA)")
            ),
            div(class = "p-6", reactableOutput("tbl_sca"))
        )
    )
  })
  
  output$tbl_f_gca <- renderReactable({ reactable(analysis_results()$f_gca, compact = TRUE) })
  output$tbl_m_gca <- renderReactable({ reactable(analysis_results()$m_gca, compact = TRUE) })
  output$tbl_sca <- renderReactable({ reactable(analysis_results()$sca, compact = TRUE) })
}

shinyApp(ui, server)
