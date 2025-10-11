# Quick Real TCGA Data Download - FIXED
# Downloads a small subset of real TCGA data

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

cat("Quick TCGA Data Download (Fixed)\n")
cat("=================================\n\n")

download_subset <- function(cancer_code, cancer_name, max_samples = 50) {
  cat(paste("\n--- Downloading", cancer_name, "---\n"))

  tryCatch({
    query <- GDCquery(
      project = paste0("TCGA-", cancer_code),
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = "Primary Tumor"
    )

    sample_list <- getResults(query)
    if (nrow(sample_list) > max_samples) {
      set.seed(42)
      selected_samples <- sample(sample_list$cases, max_samples)

      query <- GDCquery(
        project = paste0("TCGA-", cancer_code),
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        sample.type = "Primary Tumor",
        barcode = selected_samples
      )
    }

    cat(paste("  Downloading", max_samples, "samples...\n"))
    GDCdownload(query, method = "api", files.per.chunk = 5)

    cat("  Processing data...\n")
    data <- GDCprepare(query)

    # Extract expression
    counts <- assay(data, "unstranded")
    expr_matrix <- log2(counts + 1)

    # Get gene info
    gene_info <- rowData(data) %>% as.data.frame()
    gene_map <- setNames(gene_info$gene_name, rownames(gene_info))

    # Filter valid genes
    valid_genes <- !is.na(gene_map) & gene_map != ""
    expr_matrix <- expr_matrix[valid_genes, ]
    gene_symbols <- gene_map[valid_genes]

    # Handle duplicates
    dup_genes <- unique(gene_symbols[duplicated(gene_symbols)])
    keep_rows <- rep(TRUE, nrow(expr_matrix))

    for (dup_gene in dup_genes) {
      dup_indices <- which(gene_symbols == dup_gene)
      mean_expr <- rowMeans(expr_matrix[dup_indices, ], na.rm = TRUE)
      keep_idx <- dup_indices[which.max(mean_expr)]
      keep_rows[dup_indices] <- FALSE
      keep_rows[keep_idx] <- TRUE
    }

    expr_matrix <- expr_matrix[keep_rows, ]
    gene_symbols <- gene_symbols[keep_rows]
    rownames(expr_matrix) <- gene_symbols

    # Extract clinical data - FIXED to handle missing fields
    clinical <- colData(data) %>% as.data.frame()
    n_samples <- ncol(expr_matrix)

    # Create clinical data frame with safe access
    clinical_clean <- data.frame(
      sample_id = colnames(expr_matrix),
      stringsAsFactors = FALSE
    )

    # Add fields if they exist, otherwise use NA
    clinical_clean$age <- if("age_at_diagnosis" %in% names(clinical)) clinical$age_at_diagnosis else rep(NA, n_samples)
    clinical_clean$gender <- if("gender" %in% names(clinical)) clinical$gender else rep(NA, n_samples)
    clinical_clean$vital_status <- if("vital_status" %in% names(clinical)) tolower(clinical$vital_status) else rep("alive", n_samples)
    clinical_clean$days_to_death <- if("days_to_death" %in% names(clinical)) clinical$days_to_death else rep(NA, n_samples)
    clinical_clean$days_to_last_followup <- if("days_to_last_follow_up" %in% names(clinical)) clinical$days_to_last_follow_up else rep(1000, n_samples)
    clinical_clean$tumor_stage <- if("tumor_stage" %in% names(clinical)) clinical$tumor_stage else rep(NA, n_samples)

    # Create survival columns
    clinical_clean$survival_time <- ifelse(
      !is.na(clinical_clean$days_to_death),
      clinical_clean$days_to_death,
      clinical_clean$days_to_last_followup
    )

    clinical_clean$event <- ifelse(
      !is.na(clinical_clean$days_to_death) & clinical_clean$vital_status == "dead",
      1,
      0
    )

    # Filter for common cancer genes
    common_genes <- c(
      "TP53", "KRAS", "EGFR", "BRCA1", "BRCA2", "PTEN", "PIK3CA",
      "APC", "MYC", "ERBB2", "BRAF", "CDKN2A", "RB1", "NRAS",
      "FGFR1", "FGFR2", "CDK4", "MDM2", "CCND1", "VEGFA",
      "CTNNB1", "SMAD4", "ATM", "STK11", "ALK"
    )

    available_genes <- intersect(common_genes, rownames(expr_matrix))
    expr_subset <- expr_matrix[available_genes, , drop = FALSE]

    cat(paste("  ✓", cancer_name, ":", ncol(expr_subset), "samples,", nrow(expr_subset), "genes\n"))

    list(
      expression = expr_subset,
      clinical = clinical_clean
    )

  }, error = function(e) {
    cat(paste("  ✗ Error:", e$message, "\n"))
    return(NULL)
  })
}

# Download
tcga_data <- list()

tcga_data[["BRCA"]] <- download_subset("BRCA", "Breast Cancer", 50)
tcga_data[["LUAD"]] <- download_subset("LUAD", "Lung Adenocarcinoma", 50)
tcga_data[["COAD"]] <- download_subset("COAD", "Colon Adenocarcinoma", 50)

# Remove NULLs
tcga_data <- tcga_data[!sapply(tcga_data, is.null)]

# Save
if (length(tcga_data) > 0) {
  saveRDS(tcga_data, "demo_data.rds")
  cat("\n✓ Real TCGA data saved to demo_data.rds\n")
  cat("Cancer types:", paste(names(tcga_data), collapse = ", "), "\n")
  cat("Total samples:", sum(sapply(tcga_data, function(x) ncol(x$expression))), "\n")
  cat("\nYou can now run the app!\n")
  cat("Run: library(shiny); runApp()\n")
} else {
  cat("\n✗ All downloads failed\n")
}
