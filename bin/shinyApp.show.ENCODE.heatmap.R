library(shiny)
library(ggplot2)
library(reshape2)

# Input correlation matrix file path
d_mk_cor_file = '/Users/guanjuexiang/Documents/projects/git/JMnorm/docs/K562.raw_sigmat.mk.cor.txt'

# read correlation matrix
d_mk_cor = as.matrix(read.table(d_mk_cor_file, header = T, row.names = 1))

# read correlation matrix
d_mk_cor_shiny = as.matrix(d_mk_cor)
colnames(d_mk_cor_shiny) = rownames(d_mk_cor_shiny) = colnames(d_mk_cor_shiny)

# Define UI
ui <- fluidPage(
  titlePanel("ENCODE Cross-Feature Correlation Heatmap (K562)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("columns",
                  "Select Columns:",
                  choices = colnames(d_mk_cor_shiny),
                  selected = c('ATAC.seq', 'DNase-seq','H3K27ac.human','H3K27me3.human','H3K4me1.human','H3K4me3.human','H3K9me3.human','H3K36me3.human'),
                  multiple = TRUE)
    ),
    mainPanel(
      plotOutput("heatmap", width = "100%", height = "1000px")
    )
  )
)

# Define server
server <- function(input, output) {
  output$heatmap <- renderPlot({
    req(input$columns)
    
    # Subset correlation matrix
    selected_columns <- intersect(colnames(d_mk_cor_shiny), input$columns)
    subset_matrix <- d_mk_cor_shiny[selected_columns, selected_columns]
    
    # Convert matrix to long format for ggplot2
    long_data <- melt(subset_matrix)
    print(long_data)
    
    # Create heatmap
    p <- ggplot(long_data, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18), # Increase x-axis text size
        axis.text.y = element_text(size = 18), # Increase y-axis text size
        axis.title = element_text(size = 18)) + # Increase axis title size
        coord_fixed(ratio = 1) 
    
    return(p)
  })
}

# Run the application
shinyApp(ui = ui, server = server)




