library(shiny)
library(shinydashboard)
library(forestplot)
library(ggplot2)
library(dplyr)
library(DT)
library(forcats)
library(rsconnect)

#---Load data---####

#Load meta-analysis results
filenames <- list.files(path = "www/metaanalysis", pattern = "*.txt", full.names = TRUE)
data <- lapply(filenames, read.table, header = T)

#Load separate cohort results
filenames_cohorts <-list.files(path = "www/cohorts", pattern = "*.txt", full.names = TRUE)
cohorts_data <- lapply(filenames_cohorts, read.table, header = T)

#Load progression results 
aao <- read.table("www/progression/AAO_results.tab", header = T, sep = "\t")

#Need to recode the p-value threshold colunn a bit
aao <- aao %>% 
  mutate(PRS_threshold = ifelse(PRS_threshold == "Pt_5e-08", "5e-08",
                                ifelse(PRS_threshold == "Pt_9.9e-06", "1e-05",
                                       ifelse(PRS_threshold == "Pt_0.05", "0.05", NA))))

hy3 <- read.table("www/progression/HY3_results.tab", header = T, sep = "\t")

#Change the p-value threshold column to text not numeric otherwise this creates problems for the summary stats table later
hy3 <- hy3 %>% 
  mutate(PRS_threshold = as.character(PRS_threshold))

#Load separate cohort results for progression PRS
filenames_cohorts_progHY3 <- list.files(path = "www/progression", pattern = "*.txt", full.names = TRUE)
cohorts_hy3 <- lapply(filenames_cohorts_progHY3, read.table, header = T)

#Load functions
source("functions.R")

#---User interface---####

ui <- dashboardPage(
  
  #Header
  dashboardHeader(title = "Pathway-specific PRS vs clinical outcomes",
                  titleWidth = 400),
  
  #Sidebar menu
  dashboardSidebar(
    width = 400,
    sidebarMenu(
    menuItem("Pathway-PRS", tabName = "PRSwithPDloci", icon = icon("chart-bar"), selected = TRUE),
    menuItem("Progression PRS", tabName = "progressionPRS", icon = icon("chart-line")),
    menuItem("About", tabName = "about", icon = icon("question"))
    )
  ),
  
  #Main body
  dashboardBody(
    
    #Text size
    tags$style(type = "text/css",
               "label { font-size: 20px; }"
    ),
    
    tabItems(
      
      #First tab content
      tabItem(tabName = "PRSwithPDloci",
              
              #First row
              fluidRow(
                
                #First box in top row - selection of clinical outcome
                box(
                  title = div("Clinical outcomes", style = "font-size:28px;"), solidHeader = TRUE,
                  status = "warning",
                  width = 5,
                  height = 750,
                  radioButtons("outcome", 
                              label = "Choose a clinical outcome to display results for",
                              choices = c("Age at onset", 
                                          "Hoehn and Yahr 3+",
                                          "MDS-UPDRS2", 
                                          "MDS-UPDRS3",
                                          "MoCA",
                                          "Dementia",
                                          "PDQ8",
                                          "RBD"),
                              selected = "Age at onset"
                              )
                  ),
                
                #Second box in top row - display forest plot for selected outcome
                box(
                  title = textOutput("value"), solidHeader = TRUE,
                  status = "primary",
                  width = 7,
                  height = 750,
                  #Change text size of box header
                  tags$head(tags$style("#value{font-size: 28px;
                                 }")
                  ),
                  
                  #Dropdown menu for pathways
                  selectInput("pathwaySelect", label = "Choose a pathway to show the meta-analysis results with individual cohorts", 
                              choices = list("All pathways" = "all", 
                                             "Adaptive immune" = "adaptive_immune", 
                                             "Alpha synuclein" = "alpha_synuclein",
                                             "Innate immune" = "innate_immune",
                                             "Lysosome" = "lysosome",
                                             "Endocytosis" = "endocytosis",
                                             "Microglia" = "microglia",
                                             "Monocytes" = "monocytes",
                                             "Mitochondria" = "mitochondria"), 
                              selected = "All pathways"),
                  
                  #Text for selected pathway and selected outcome
                  column(12, align = "center", htmlOutput("plotTitleHTML"), 
                         #Change text size of box header
                         tags$head(tags$style("#plotTitleHTML{font-size: 22px;
                                 }")),
                         ), 
                  
                  #Forest plot for selected pathway, or all pathways
                  plotOutput("plot", height = 540)
                )
              ),
              
              #Second row
              fluidRow(
                
                #Second row, single box with full page width - summary statistics
                box(width = 12, title = "Summary statistics",
                    DTOutput("table")
                    )
              )
      ),
      
      #Second tab content
      tabItem(tabName = "progressionPRS",
              
              #First row
              fluidRow(
                
                #First box in top row - selection of clinical outcome
                box(
                  title = div("Clinical outcomes", style = "font-size:28px;"), solidHeader = TRUE,
                  status = "warning",
                  width = 5,
                  height = 750,
                  radioButtons("progression_outcome", 
                               label = "Choose a clinical outcome to display results for",
                               choices = c("Age at onset", 
                                           "Hoehn and Yahr 3+"),
                               selected = "Age at onset"
                  )
                ),
                
                #Second box in top row - display forest plot for selected outcome
                box(
                  title = textOutput("progression_value"), solidHeader = TRUE,
                  status = "primary",
                  width = 7,
                  height = 750,
                  #Change text size of box header
                  tags$head(tags$style("#progression_value{font-size: 28px;
                                 }")
                  ),
                  
                  #Dropdown menu for pathways
                  selectInput("progression_pval", label = "Choose a PRS p-value threshold (pathway-specific PRSs were not analysed for progression PRSs)", 
                              choices = list(
                                "p < 5E-08" = "5e-08",
                                "p < 1E-05" = "1e-05",
                                "p < 0.05" = "0.05"), 
                              selected = "p < 5E-08"
                  ),
                  
                  #Text for selected p-value threshold and selected clinical outcome
                  column(12, align = "center", htmlOutput("progression_plotTitleHTML"), 
                         #Change text size of box header
                         tags$head(tags$style("#progression_plotTitleHTML{font-size: 28px;
                                 }")),
                  ), 
                  
                  #Forest plot for selected  PRS p-value threshold
                  plotOutput("progression_plot", height = 540)
                )
              ),
              
              #Second row
              fluidRow(
                
                #Second row, single box with full page width - summary statistics
                box(width = 12, title = "Summary statistics",
                    DTOutput("progression_table")
                )
              )
              
      ),
      
      
      #Third tab content
      tabItem(tabName = "about",
              title = div("About", style = "font-size:28px;"), solidHeader = TRUE,
              
              #Text
              div(
                style = "font-weight:bold;
                     font-size:130%;
                    margin-bottom:12px",
                "For more information on this data, please see",
                a(
                  "Tan et al. 2024 (under review)",
                  href = "https://www.biorxiv.org/content/10.1101/585836v2",
                  target = "_blank"
                ),
                "or contact:"
              ),
              div(span("Manuela Tan", style = "font-weight:bold; font-size:110%"),
                  div(a(href="mailto:manuela.tan@medisin.uio.no", "manuela.tan@medisin.uio.no")),
                  div("Department of Neurology"),
                  div("Oslo University Hospital"),
                  div("Oslo"),
                  div("Norway"))
      )
    )
  )
)



#---Server logic---####

server <- function(input, output) {
  
  ## FOR PATHWAY-SPECIFIC PRS (FIRST TAB) ##
  
  #Select relevant outcome dataset - all pathways vs. selected outcome
  dataInput <- reactive({
    
    #Switch text values for numeric values
    #These are in the order that the datasets are read in (alphabetical order)
    #The data is a list of dataframes, just need to index the correct dataframe
    data[[switch(input$outcome,
                   "Age at onset" = 1, 
                   "Dementia" = 2,
                   "Hoehn and Yahr 3+" = 3,
                   "MDS-UPDRS2" = 4,
                   "MDS-UPDRS3" = 5,
                   "MoCA" = 6,
                   "PDQ8" = 7,
                   "RBD" = 8)]]
  })
  
  #Select relevant pathway dataset - this is the meta-analysis result for the specific pathway selected
  meta_pathway <- reactive({
    
    #If pathway selected is 'all' just return original dataInput with all pathway results
    if(input$pathwaySelect == "all"){
      return(dataInput())
      
    } else {
      
      #If a pathway is selected, filter dataframe for selected pathway
      meta <- dataInput()
      meta <- meta %>%
        filter(pathway == input$pathwaySelect) %>%
        mutate(dataset = "Random effects meta-analysis",
               N = NA) %>%
        select(dataset, pathway, Coeff = RE_SMD, se = RE_se, Pvalue = pval, N)
    }
  })
  
  
  #Select relevant cohort dataset - selected outcome
  cohort_outcome <- reactive({
    
    #Index the cohorts data for the correct outcome
    cohorts_data[[switch(input$outcome,
                           "Age at onset" = 1, 
                           "Dementia" = 2,
                           "Hoehn and Yahr 3+" = 3,
                           "MDS-UPDRS2" = 4,
                           "MDS-UPDRS3" = 5,
                           "MoCA" = 6,
                           "PDQ8" = 7,
                           "RBD" = 8)]]
  })

  #Select relevant pathway dataset - separate cohort results
  cohort_outcome_pathway <- reactive({
    
    if(input$pathwaySelect == "all"){
      NULL
      
    } else {
      
      cohorts <- cohort_outcome()
      
      #Filter for selected pathway
      data_filtered <- cohorts %>%
        filter(pathway == input$pathwaySelect)
    }
  })
  
  #For selected pathway, combine the meta-analysis result with the separate cohort results
  #So we can plot the forest plot with cohorts separated
  finalInput <- reactive({
    
    if(input$pathwaySelect == "all"){
      NULL
    } else {
      meta_results <- meta_pathway()
      cohort_results <- cohort_outcome_pathway()
      rbind(cohort_results, meta_results)
    }
    
    
  })
  
  
  ## FOR PROGRESSION PRS (SECOND TAB) ##
  
  #Select relevant progression PRS dataset dependent on p-value threshold selected
  progCohorts_Input <- reactive({
    
    #Switch is used to swap text values for numeric values
    #The data is a list of dataframes, just need to index the correct dataframe
    cohorts <- cohorts_hy3[[switch(input$progression_pval,
                                   "0.05" = 1,
                                   "1e-05" = 2,
                                   "5e-08" = 3)]]

    return(cohorts)
  })
  
  
  #Select the relevant progression PRS meta-analysis results dependent on the outcome selected
  progMeta_Input <- reactive({
    
    #Filter progression PRS meta-analysis results for appropriate p-value threshold
    meta <- if (input$progression_outcome == "Age at onset") {
      aao %>% filter(PRS_threshold == input$progression_pval)
    } else {
      hy3 %>% filter(PRS_threshold == input$progression_pval)
    }
    return(meta)
  })
  
  
  output$debug_progCohorts_Input <- renderPrint({
    progCohorts_Input()
  })
  
  output$debug_progMeta_Input <- renderPrint({
    progMeta_Input()
  })
  
  #Final progression dataset
  progFinal_Input <- reactive({
    final <- if (input$progression_outcome == "Age at onset") {
      
      #If the outcome+PRS is age at onset, there is no meta-analysis so just return and tidy the dataset filtered by p-value threshold
      foxinsight_forplot <- progMeta_Input()
      
      #Tidy the dataframe a bit - this is for the forest plot
      foxinsight_forplot <- foxinsight_forplot %>% 
        mutate(dataset = "Fox Insight") %>% 
        rename(Pvalue = p.value) %>%
        select(dataset, Coeff = Beta, se = SE, Pvalue, N = N_inds, pval_threshold = PRS_threshold)
      
    } else {
      
      meta_results <- progMeta_Input()
      cohort_results <- progCohorts_Input()

      # Rename and select necessary columns
      meta_results <- meta_results %>%
          mutate(dataset = "Random effects meta-analysis") %>% 
          rename(Pvalue = p.value) %>%
          select(dataset, Coeff = Beta, se = SE, Pvalue, N = N_inds, pval_threshold = PRS_threshold)

      # Combine individual cohort summary statistics and meta-analysis results for the forest plot
      combined_results <- rbind(cohort_results, meta_results)
    }
    
    return(final)
  })

    
  #Title for forest plot according to selected outcome
  output$value <- renderText({
    paste("Forest plot for pathway-specific PRSs vs.", input$outcome)
  })
  
  #Text header for forest plot according to selected outcome and pathway
  #Switch back to original values so can display full text
  output$plotTitleHTML <- renderText({
    paste("<b>",switch(input$pathwaySelect,
                 "all" = "All pathways", 
                 "adaptive_immune" = "Adaptive immune", 
                 "alpha_synuclein" = "Alpha synuclein",
                 "innate_immune" = "Innate immune",
                 "lysosome" = "Lysosome",
                 "endocytosis" = "Endocytosis",
                 "microglia" = "Microglia",
                 "monocytes" = "Monocytes",
                 "mitochondria" = "Mitochondria"), "vs", input$outcome, "</b>", sep = " ")
  })
  
  #Make forest plot for pathway PRSs
  output$plot <- renderPlot({
    
    #makeForestPlot and makeForestPlot_cohorts are separate custom functions stored in functions.R
    #If pathways = 'all', show the forest plot with meta-analysis results for all the pathways
    if(input$pathwaySelect == "all"){
      makeForestPlot(dataInput())
    } else {
      #If a specific pathway is selected, plot the forest plot with separate cohort results for selected outcome and selected pathway
      makeForestPlot_cohorts(finalInput())
    }
    
  })
  
  
  #Make summary stats table using DT package
  output$table <- renderDT({
    table <- dataInput()
    #Assign column names to table
    colnames(table) <- c("Pathway", "Beta", "SE", "95% CI lower", "95% CI upper", "z", "p-value", "N_studies", "I2", "Cochrans Q p-val", "N_inds")
    
  
    #Round numeric columns to 3 digits
    table <- table %>% 
      mutate_if(is.numeric, round, digits = 3)
    
    #Output table
    datatable(table, rownames = FALSE, options = list(
              
              #Other options to display download buttons
              paging = TRUE,
              searching = TRUE,
              fixedColumns = TRUE,
              autoWidth = TRUE,
              ordering = TRUE,
              dom = 'tB',
              buttons = c('csv', 'excel')),
    
    extensions = 'Buttons',
    class = "display")
  })
  
  
  #Make forest plot for progression PRSs
  output$progression_plot <- renderPlot({
    
    #If the outcome/PRS is age at onset, just make the forest plot without meta-analysis results as this is just Fox Insight dataset
    if(input$progression_outcome == "Age at onset"){
      makeProgressionPRSplot_nometa(progFinal_Input())
    } else {
      #If other outcome/PRS is selected, make the forest plot with each cohort and meta-analysis results
      makeProgressionPRSplot(progFinal_Input())
    }
    
  })
  
  output$progression_table <- renderDT({
    table <- progMeta_Input()
    #Assign column names to table
    colnames(table) <- c("PRS p-value threshold", "Beta", "SE", "95% CI lower", "95% CI upper", "z", "p-value", "N_studies", "I2", "Cochrans Q p-val", "N_inds")
    
    #Round numeric columns to 3 digits
    table <- table %>% 
      mutate_if(is.numeric, round, digits = 3)
    
    #Code to explicitly show NAs in the table (instead of just blank/empty cells)
    #This shows the NAs as italic and grey text
    rowCallback <- c(
      "function(row, data){",
      "  for(var i=0; i<data.length; i++){",
      "    if(data[i] === null){",
      "      $('td:eq('+i+')', row).html('NA')",
      "        .css({'color': 'rgb(151,151,151)', 'font-style': 'italic'});",
      "    }",
      "  }",
      "}"  
    )
    
    #Output table
    datatable(table, rownames = FALSE, options = list(
      rowCallback = JS(rowCallback), #This is to display the NAs as text
      #Other options to display download buttons
      paging = TRUE,
      searching = TRUE,
      fixedColumns = TRUE,
      autoWidth = TRUE,
      ordering = TRUE,
      dom = 'tB',
      buttons = c('csv', 'excel')),
      
      extensions = 'Buttons',
      class = "display")
              
              
  })
  
  
}





#---Run the app---####

shinyApp(ui, server)