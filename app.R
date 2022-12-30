library(shiny)
library(shinyMobile)
library(dplyr)
library(ggplot2)
library(tidyr)
library(DT)

data_page <- tabPanel(
  title = 'Data input',
  sidebarLayout(
    sidebarPanel(
      # Select csv file
      fileInput(inputId = 'raw_file', label = 'Select .csv file',
                multiple = FALSE,
                accept = c('.csv')),
      

    ),
    mainPanel(
      # Define treatments
      selectInput(inputId = 'control', label ='Control:',
                  choices = NULL,
                  multiple = TRUE),
      
      selectInput(inputId = 'treatment', label ='Treatment:',
                  choices = NULL,
                  multiple = TRUE)
    )
  )
)

processing_page <- tabPanel(
  title = 'Data Processing',
  sidebarLayout(
    
    sidebarPanel(
      
      # Inputs
      #Select filtering threshhold
      sliderInput(inputId = 'filter', label = '% of identification',
                  min = 0, max = 100,
                  value = 0),
      
      # Select missing value imputation
      selectInput(inputId = 'imputation', label = 'MV imputation',
                  choices = c('none'= 'none',
                              '0.001' = '0.001')),
      
      # Select transformation
      selectInput(inputId = 'transformation', label = 'Transformation',
                  choices = c('none' = 'none',
                              'log2' = 'log2',
                              'log10'= 'log10')),
      
      # Select normalization
      selectInput(inputId = 'normalization', label = 'Normalization',
                  choices = c('none' = 'none',
                              'z-score' = 'z-score')),

      downloadButton(outputId = 'download', label = 'Download Processed Data')
      
    ),
    
    # Main panel for displaying outputs
    mainPanel(
      
      tabsetPanel(type = "tabs",
                  
                  tabPanel('Processed Data', dataTableOutput('processed_data')),
                  
                  tabPanel('Missing Data',fluidPage(
                    tags$style(type = "text/css", "#missing_map {height: 80vh !important;}"),
                    plotOutput('missing_map'))),
            
                  tabPanel('Visual inspection', plotOutput('norm_plot')))
      
    )
  )
)

plots_page <- tabPanel(
  title = 'Plots',
  tabsetPanel(type = 'tabs',
              tabPanel('Boxplots', fillPage(
                tags$style(type = "text/css", "#boxplots {height: 80vh !important;}"),
                plotOutput('boxplots'))),
              tabPanel('Heatmap', 
                       sidebarPanel(
                         ## PUT CLUSTERS and METHOD ##
                       ),
                       mainPanel(
                         fillPage(
                           tags$style(type = "text/css", "#heatmap {height: 80vh !important;}"),
                           plotOutput('heatmap')))
              ))
                       
)


ui <- navbarPage(
  title = 'Omics Analyzer',
  data_page,
  processing_page,
  plots_page
)

# Functions

# Data processes

# Data filtering
filtering_function <- function(df, filter_value){
  
  filter_value <- 1 - filter_value/100
  
  # Calculate the percentage of zero values for each group
  zero_percentages <- df %>%
    group_by(ID) %>% 
    summarize(zero_percentage = sum(INTENSITY == 0) / n())
  
  # Filter the data frame to include only the rows where the percentage of zero values is greater than or equal to the specified percent
  df_filtered <- left_join(df, zero_percentages, by = c("ID")) %>%
    filter(zero_percentage <= filter_value) %>%
    select(-zero_percentage)
    
  return(df_filtered)
}

imputation_function <- function(df_filtered, impute_value){
  if (impute_value == 'none'){
    df_imputed <- df_filtered
  } else if (impute_value == '0.001'){
    df_imputed <- df_filtered %>%
      mutate(INTENSITY = ifelse(INTENSITY == 0, 0.001, INTENSITY))
  }
  
  return(df_imputed)
}

# Data transformation
transform_function <- function(df_imputed, transform_value){
  if (transform_value == 'none'){
    df_transform <- df_imputed
  } else if (transform_value == 'log2'){
    df_transform <- df_imputed %>%
      mutate(INTENSITY = log2(INTENSITY))
  } else if (transform_value == 'log10'){
    df_transform <- df_imputed %>%
      mutate(INTENSITY = log(INTENSITY))
  }
  
  return(df_transform)
}
  
# Data normalization
normalize_function <- function(df_transform, norm_value){
  if (norm_value == 'none'){
    df_norm <- df_transform
  } else if (norm_value == 'z-score'){
    df_norm <- df_transform %>%
      group_by(ID, TREATMENT) %>%
      mutate(INTENSITY = scale(INTENSITY))
  }
  
  return(df_norm)
}

wide_to_long_function <- function(df_wide, control_vector, treatment_vector){
  temp <-tibble::rownames_to_column(df_wide, 'ID')
  df_long <- temp %>%
    pivot_longer(cols = !c(1), names_to = 'SAMPLE', values_to = 'INTENSITY') %>%
    mutate(TREATMENT = case_when(SAMPLE %in% control_vector ~ 'Control',
                                 SAMPLE %in% treatment_vector ~ 'Treatment'))
  return(df_long)
}

processing_function <- function(df, filter_val, impute_val, transformation_val, normalization_val) {
  df_filtered <- filtering_function(df, filter_val)
  df_imputed <- imputation_function(df_filtered, impute_val)
  df_transform <- transform_function(df_imputed, transformation_val)
  df_norm <- normalize_function(df_transform, normalization_val)
  
  return(df_norm)
}

server <- function(input, output, session){
  
  # Read input dataframe
  df <- reactive({
    req(input$raw_file)
    read.csv(input$raw_file$datapath, row.names = 1)
  })
  
  observe({
    updateSelectInput(session, 'control', choices = colnames(df()))
    updateSelectInput(session, 'treatment', choices = colnames(df())) 
  })
  
  
  # Processed Data panel
  output$processed_data <- renderDataTable({
    datatable(df_work(), options = list(scrollY = "400px"))
  })
  
  # # Reactive input variables
  df_long <- reactive(wide_to_long_function(df(), control_vector(), treatment_vector()))
  df_work <- reactive(processing_function(df_long(), filter_value(), impute_value(), transformation_value(), normalization_value()))
  
  control_vector <- reactive(input$control)
  treatment_vector <- reactive(input$treatment)
  filter_value <- reactive(input$filter)
  impute_value <- reactive(input$imputation)
  transformation_value <- reactive(input$transformation)
  normalization_value <- reactive(input$normalization)
  # df_wide <- reactive(processes_function(df(), 
  #                                         filter_value(), 
  #                                         transformation_value(), 
  #                                         normalization_value()))
  
  # Missing Data panel
  output$missing_map <- renderPlot({
    df_work() %>%
      ggplot(aes(x = SAMPLE, y = ID, fill = ifelse(INTENSITY != 0, 'Non-missing', 'Missing'))) +
      geom_tile(color = 'black') +
      scale_fill_manual(values = c("Missing" = "red", "Non-missing" = "green")) +
      theme_linedraw() +
      theme(legend.title = element_blank(), axis.text.x = element_text(angle = 90)) +
      xlab(NULL) +
      ylab(NULL) 
  })
  
  
  output$norm_plot <- renderPlot({
    df_work() %>%
      ggplot(aes(x = SAMPLE, y = INTENSITY, fill = TREATMENT))+
      geom_boxplot(position = )+
      geom_point()+
      theme_linedraw()+
      theme(legend.position = 'none', text = element_text(size=20))+
      xlab(NULL)
  })
  
  output$download <- downloadHandler( 
    filename <- function() {'processed_data.csv'},
    content <- function(file) {
      write.csv(df_work(), file, row.names = FALSE)
    } 
  )
  
  # # Plots Panels
  output$boxplots <- renderPlot({
    df_work() %>%
      ggplot(aes(x = TREATMENT, y = INTENSITY, color = TREATMENT))+
      geom_boxplot()+
      facet_wrap(~ID)+
      theme_classic()+
      theme(text = element_text(size=20), axis.text.x = element_text(angle = 90))+
      xlab(NULL)
  })

  output$heatmap <- renderPlot({
    
  })
    
}

shinyApp(ui = ui, server = server)
