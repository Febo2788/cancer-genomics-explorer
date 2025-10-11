# Cancer Genomics Explorer - Shiny App
# Interactive visualization of TCGA cancer genomics data

library(shiny)
library(ggplot2)
library(dplyr)
library(survival)
library(survminer)
library(tidyr)

# Load data
if (file.exists("demo_data.rds")) {
  tcga_data <- readRDS("demo_data.rds")
} else {
  stop("Data file not found. Please run generate_demo_data.R first.")
}

# Extract available cancer types and genes
cancer_types <- names(tcga_data)
all_genes <- rownames(tcga_data[[1]]$expression)

# UI
ui <- fluidPage(
  titlePanel("Cancer Genomics Explorer"),

  sidebarLayout(
    sidebarPanel(
      h4("Data Selection"),
      selectInput(
        "cancer_type",
        "Cancer Type:",
        choices = setNames(
          cancer_types,
          c("Breast Cancer (BRCA)",
            "Lung Adenocarcinoma (LUAD)",
            "Colon Adenocarcinoma (COAD)")
        ),
        selected = cancer_types[1]
      ),

      selectInput(
        "gene",
        "Gene:",
        choices = all_genes,
        selected = "TP53"
      ),

      hr(),

      h4("Visualization Options"),
      radioButtons(
        "plot_type",
        "Expression Plot Type:",
        choices = c("Boxplot" = "box", "Violin Plot" = "violin"),
        selected = "box"
      ),

      sliderInput(
        "survival_cutoff",
        "Survival Analysis Cutoff (percentile):",
        min = 25,
        max = 75,
        value = 50,
        step = 5
      ),

      hr(),

      h4("About"),
      p("This app visualizes cancer genomics data from TCGA."),
      p("Explore gene expression patterns, survival outcomes, and correlations."),
      tags$small("Data source: The Cancer Genome Atlas (TCGA)")
    ),

    mainPanel(
      tabsetPanel(
        id = "tabs",

        tabPanel(
          "Expression",
          h3("Gene Expression Across Cancer Types"),
          plotOutput("expression_plot", height = "400px"),
          downloadButton("download_expression", "Download Plot"),
          hr(),
          verbatimTextOutput("expression_summary")
        ),

        tabPanel(
          "Survival",
          h3("Survival Analysis"),
          p("Kaplan-Meier curves comparing high vs. low gene expression"),
          plotOutput("survival_plot", height = "500px"),
          downloadButton("download_survival", "Download Plot"),
          hr(),
          verbatimTextOutput("survival_summary")
        ),

        tabPanel(
          "Correlations",
          h3("Gene Correlations"),
          p("Top genes correlated with selected gene"),
          plotOutput("correlation_plot", height = "500px"),
          downloadButton("download_correlation", "Download Plot"),
          hr(),
          tableOutput("correlation_table")
        ),

        tabPanel(
          "Data Info",
          h3("Dataset Information"),
          verbatimTextOutput("data_info")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {

  # Reactive data based on selected cancer type
  current_data <- reactive({
    req(input$cancer_type)
    tcga_data[[input$cancer_type]]
  })

  # Expression data for selected gene
  gene_expression <- reactive({
    req(input$gene)
    data <- current_data()
    expr_values <- data$expression[input$gene, ]

    data.frame(
      sample = names(expr_values),
      expression = as.numeric(expr_values),
      cancer_type = input$cancer_type
    )
  })

  # Expression plot across all cancer types
  output$expression_plot <- renderPlot({
    req(input$gene)

    # Combine data from all cancer types
    all_expr <- lapply(cancer_types, function(ct) {
      expr_values <- tcga_data[[ct]]$expression[input$gene, ]
      data.frame(
        expression = as.numeric(expr_values),
        cancer_type = ct,
        stringsAsFactors = FALSE
      )
    })
    plot_data <- bind_rows(all_expr)

    # Create plot
    p <- ggplot(plot_data, aes(x = cancer_type, y = expression, fill = cancer_type)) +
      labs(
        title = paste("Expression of", input$gene, "Across Cancer Types"),
        x = "Cancer Type",
        y = "Expression (log2 FPKM + 1)",
        fill = "Cancer Type"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      scale_fill_brewer(palette = "Set2")

    if (input$plot_type == "box") {
      p + geom_boxplot(alpha = 0.7) +
          geom_jitter(width = 0.2, alpha = 0.3, size = 0.5)
    } else {
      p + geom_violin(alpha = 0.7) +
          geom_boxplot(width = 0.1, fill = "white", alpha = 0.8)
    }
  })

  # Expression summary statistics
  output$expression_summary <- renderText({
    data <- gene_expression()

    paste0(
      "Gene: ", input$gene, "\n",
      "Cancer Type: ", input$cancer_type, "\n",
      "Number of Samples: ", nrow(data), "\n",
      "Mean Expression: ", round(mean(data$expression, na.rm = TRUE), 3), "\n",
      "Median Expression: ", round(median(data$expression, na.rm = TRUE), 3), "\n",
      "SD: ", round(sd(data$expression, na.rm = TRUE), 3), "\n",
      "Range: ", round(min(data$expression, na.rm = TRUE), 3), " - ",
      round(max(data$expression, na.rm = TRUE), 3)
    )
  })

  # Survival analysis
  survival_data <- reactive({
    req(input$gene)

    tryCatch({
      data <- current_data()
      clinical <- data$clinical
      expr_values <- data$expression[input$gene, ]

      # Match samples
      common_samples <- intersect(names(expr_values), clinical$sample_id)

      if (length(common_samples) == 0) {
        return(NULL)
      }

      clinical <- clinical[clinical$sample_id %in% common_samples, ]
      expr_values <- expr_values[common_samples]

      # Ensure survival columns exist
      if (!"survival_time" %in% names(clinical)) {
        clinical$survival_time <- rep(1000, nrow(clinical))
      }
      if (!"event" %in% names(clinical)) {
        clinical$event <- rep(0, nrow(clinical))
      }

      # Categorize into high/low expression
      cutoff <- quantile(expr_values, probs = input$survival_cutoff / 100, na.rm = TRUE)
      clinical$expression_group <- ifelse(expr_values > cutoff, "High", "Low")
      clinical$expression_value <- expr_values

      clinical
    }, error = function(e) {
      return(NULL)
    })
  })

  output$survival_plot <- renderPlot({
    req(input$gene)

    tryCatch({
      surv_data <- survival_data()

      # Check if data is available
      if (is.null(surv_data)) {
        plot.new()
        text(0.5, 0.5, "No survival data available", cex = 1.5)
        return()
      }

      # Remove samples with missing survival data
      surv_data <- surv_data[!is.na(surv_data$survival_time) &
                              !is.na(surv_data$event), ]

      if (nrow(surv_data) < 10) {
        plot.new()
        text(0.5, 0.5, "Insufficient survival data available", cex = 1.5)
        return()
      }

      # Fit survival curves directly using formula with column names
      fit <- survfit(Surv(survival_time, event) ~ expression_group, data = surv_data)

      # Plot with survminer
      ggsurvplot(
        fit,
        data = surv_data,
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.height = 0.25,
        xlab = "Time (days)",
        ylab = "Survival Probability",
        title = paste("Survival Analysis:", input$gene, "in", input$cancer_type),
        legend.title = "Expression",
        legend.labs = c("High", "Low"),
        palette = c("#E7B800", "#2E9FDF"),
        ggtheme = theme_minimal()
      )
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Error:", e$message), cex = 1.2)
    })
  })

  # Survival summary
  output$survival_summary <- renderText({
    surv_data <- survival_data()

    if (is.null(surv_data)) {
      return("No survival data available")
    }

    surv_data <- surv_data[!is.na(surv_data$survival_time) &
                            !is.na(surv_data$event), ]

    if (nrow(surv_data) < 10) {
      return("Insufficient data for survival analysis")
    }

    high_group <- surv_data[surv_data$expression_group == "High", ]
    low_group <- surv_data[surv_data$expression_group == "Low", ]

    paste0(
      "Survival Analysis Summary\n",
      "-------------------------\n",
      "Cutoff: ", input$survival_cutoff, "th percentile\n\n",
      "High Expression Group:\n",
      "  Samples: ", nrow(high_group), "\n",
      "  Deaths: ", sum(high_group$event), "\n",
      "  Median survival: ", round(median(high_group$survival_time, na.rm = TRUE)), " days\n\n",
      "Low Expression Group:\n",
      "  Samples: ", nrow(low_group), "\n",
      "  Deaths: ", sum(low_group$event), "\n",
      "  Median survival: ", round(median(low_group$survival_time, na.rm = TRUE)), " days"
    )
  })

  # Correlation analysis
  correlation_data <- reactive({
    req(input$gene)
    data <- current_data()
    expr_matrix <- data$expression

    # Get expression of selected gene
    target_gene_expr <- expr_matrix[input$gene, ]

    # Calculate correlations with all genes
    correlations <- apply(expr_matrix, 1, function(gene_expr) {
      cor(target_gene_expr, gene_expr, use = "complete.obs")
    })

    # Create data frame
    cor_df <- data.frame(
      gene = names(correlations),
      correlation = correlations,
      stringsAsFactors = FALSE
    )

    # Remove the query gene itself
    cor_df <- cor_df[cor_df$gene != input$gene, ]

    # Sort by absolute correlation
    cor_df <- cor_df[order(abs(cor_df$correlation), decreasing = TRUE), ]

    cor_df
  })

  output$correlation_plot <- renderPlot({
    cor_data <- correlation_data()

    # Top 20 correlations
    top_cors <- head(cor_data, 20)

    # Plot
    ggplot(top_cors, aes(x = reorder(gene, correlation), y = correlation,
                         fill = correlation > 0)) +
      geom_col(alpha = 0.8) +
      coord_flip() +
      labs(
        title = paste("Top 20 Genes Correlated with", input$gene),
        subtitle = paste("in", input$cancer_type),
        x = "Gene",
        y = "Pearson Correlation Coefficient",
        fill = "Direction"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom"
      ) +
      scale_fill_manual(
        values = c("TRUE" = "#00BA38", "FALSE" = "#F8766D"),
        labels = c("Negative", "Positive")
      )
  })

  output$correlation_table <- renderTable({
    cor_data <- correlation_data()
    top_10 <- head(cor_data, 10)
    top_10$correlation <- round(top_10$correlation, 3)
    top_10
  }, striped = TRUE, hover = TRUE)

  # Dataset information
  output$data_info <- renderText({
    info <- paste0(
      "Cancer Genomics Explorer\n",
      "========================\n\n",
      "Available Cancer Types:\n"
    )

    for (ct in cancer_types) {
      n_samples <- ncol(tcga_data[[ct]]$expression)
      n_genes <- nrow(tcga_data[[ct]]$expression)
      info <- paste0(
        info,
        "  - ", ct, ": ", n_samples, " samples, ", n_genes, " genes\n"
      )
    }

    info <- paste0(
      info,
      "\nTotal Genes: ", length(all_genes), "\n",
      "\nData Source: The Cancer Genome Atlas (TCGA)\n",
      "Expression values: log2(FPKM + 1)\n",
      "\nFeatures:\n",
      "  - Gene expression visualization across cancer types\n",
      "  - Kaplan-Meier survival analysis\n",
      "  - Gene correlation analysis\n"
    )

    info
  })

  # Download handlers
  output$download_expression <- downloadHandler(
    filename = function() {
      paste0(input$gene, "_expression_", Sys.Date(), ".png")
    },
    content = function(file) {
      # Recreate the plot
      all_expr <- lapply(cancer_types, function(ct) {
        expr_values <- tcga_data[[ct]]$expression[input$gene, ]
        data.frame(
          expression = as.numeric(expr_values),
          cancer_type = ct,
          stringsAsFactors = FALSE
        )
      })
      plot_data <- bind_rows(all_expr)

      p <- ggplot(plot_data, aes(x = cancer_type, y = expression, fill = cancer_type)) +
        labs(
          title = paste("Expression of", input$gene, "Across Cancer Types"),
          x = "Cancer Type",
          y = "Expression (log2 FPKM + 1)",
          fill = "Cancer Type"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold")
        ) +
        scale_fill_brewer(palette = "Set2")

      if (input$plot_type == "box") {
        p <- p + geom_boxplot(alpha = 0.7) + geom_jitter(width = 0.2, alpha = 0.3, size = 0.5)
      } else {
        p <- p + geom_violin(alpha = 0.7) + geom_boxplot(width = 0.1, fill = "white", alpha = 0.8)
      }

      ggsave(file, plot = p, width = 10, height = 6, dpi = 300)
    }
  )

  output$download_survival <- downloadHandler(
    filename = function() {
      paste0(input$gene, "_survival_", Sys.Date(), ".png")
    },
    content = function(file) {
      tryCatch({
        surv_data <- survival_data()
        surv_data <- surv_data[!is.na(surv_data$survival_time) &
                                !is.na(surv_data$event), ]

        if (nrow(surv_data) < 10) {
          stop("Insufficient survival data")
        }

        fit <- survfit(Surv(survival_time, event) ~ expression_group, data = surv_data)

        p <- ggsurvplot(
          fit,
          data = surv_data,
          pval = TRUE,
          conf.int = TRUE,
          risk.table = TRUE,
          risk.table.height = 0.25,
          xlab = "Time (days)",
          ylab = "Survival Probability",
          title = paste("Survival Analysis:", input$gene, "in", input$cancer_type),
          legend.title = "Expression",
          legend.labs = c("High", "Low"),
          palette = c("#E7B800", "#2E9FDF"),
          ggtheme = theme_minimal()
        )

        ggsave(file, print(p), width = 10, height = 8, dpi = 300)
      }, error = function(e) {
        # Create error message plot
        png(file, width = 800, height = 600)
        plot.new()
        text(0.5, 0.5, paste("Error generating plot:", e$message))
        dev.off()
      })
    }
  )

  output$download_correlation <- downloadHandler(
    filename = function() {
      paste0(input$gene, "_correlations_", Sys.Date(), ".png")
    },
    content = function(file) {
      cor_data <- correlation_data()
      top_cors <- head(cor_data, 20)

      p <- ggplot(top_cors, aes(x = reorder(gene, correlation), y = correlation,
                                 fill = correlation > 0)) +
        geom_col(alpha = 0.8) +
        coord_flip() +
        labs(
          title = paste("Top 20 Genes Correlated with", input$gene),
          subtitle = paste("in", input$cancer_type),
          x = "Gene",
          y = "Pearson Correlation Coefficient",
          fill = "Direction"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          legend.position = "bottom"
        ) +
        scale_fill_manual(
          values = c("TRUE" = "#00BA38", "FALSE" = "#F8766D"),
          labels = c("Negative", "Positive")
        )

      ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)
