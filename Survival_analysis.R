# setwd("~/Documents/TCGA-BRCA")
#script to perform survival analysis from cbiportal TCGA-BRCA data

library(readr)
library(dplyr)

#-------------STEP:1 DATA DOWNLOAD----------------

#load gene expression data and clinical data
gene_data <- read_tsv("gene_expression.txt")
clin_data <- read_tsv("clinical_data.tsv")

colnames(clin_data)

#-----------------STEP:2 DATA PREPARARTION---------------------

# Remove rows with Nan values from gene expression data
gene <- gene_data %>% na.omit

#Extract the survival data (time to event) and event status from the clinical file. 
any(colnames(clin_data) %in% c("Sample ID", "Overall Survival (Months)","Overall Survival Status"))
which(colnames(clin_data) %in% c("Sample ID","Overall Survival (Months)","Overall Survival Status"))
clin_data[,c(3,35,36)]

# Extract survival data
survival_data <- clin_data %>%
  select("Sample ID","Overall Survival (Months)","Overall Survival Status")

# Rename columns for clarity
colnames(survival_data) <- c("Sample_ID","OS_Months", "OS_status")

# Preprocess survival status to binary (0 for living, 1 for deceased)
survival_data <-survival_data %>%
  mutate(OS_status = ifelse(OS_status == "0:LIVING", 0, 1))

str(gene)
str(survival_data)

#---------------STEP:3 FEATURE SELECTION--------------------------

# Transpose gene expression data so that sample IDs are in a single column
library(tidyr)
gene_transposed <- gene %>%
  gather(key = "Sample_ID", value = "Expression", -Hugo_Symbol, -Entrez_Gene_Id)

# Merge survival data and gene_transposed based on Sample_ID
merged_data <- survival_data %>%
  inner_join(gene_transposed, by = "Sample_ID")

colnames(merged_data)
str(merged_data)
summary(merged_data)

library(survival)

# Ensure OS_Months is numeric
merged_data$OS_Months <- as.numeric(merged_data$OS_Months)

# Check the preprocessed data
str(merged_data)

# Function to fit Cox model for a given gene and return p-value
fit_cox_model <- function(gene) {
  cox_model <- coxph(Surv(OS_Months, OS_status) ~ Expression, data = gene)
  summary(cox_model)$coefficients[5]
}

# Split data by gene
gene_list <- split(merged_data, merged_data$Hugo_Symbol)

# Fit Cox model and extract p-values for each gene
p_values <- sapply(gene_list, fit_cox_model)

# Create a data frame with genes and their p-values
gene_p_values <- data.frame(
  Hugo_Symbol = names(p_values),
  p_value = p_values
)

# Check the top genes with the lowest p-values
head(gene_p_values)

# Select top 100 genes based on p-value
top_genes <- gene_p_values %>%
  arrange(p_value) %>%
  head(100)

# Check the top 100 genes
print(top_genes)

# Load library for survival curves
library(survminer)

# Function to plot survival curves for a given gene
plot_survival_curve <- function(gene_name, data) {
  gene_data <- data %>% filter(Hugo_Symbol == gene_name)
  fit <- survfit(Surv(OS_Months, OS_status) ~ Expression > median(Expression), data = gene_data)
  ggsurvplot(fit, data = gene_data, risk.table = TRUE, pval = TRUE, 
             title = paste("Survival curve for", gene_name))
}

# Plot survival curves for the top 5 genes
for (i in 1:5) {
  plot_survival_curve(top_genes$Hugo_Symbol[i], merged_data)
}


#---------------------STEP 4: CLUSTERING--------------------------
# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(factoextra)

# Filter the merged data to include only the top 100 genes
top_gene_names <- top_genes$Hugo_Symbol
top_gene_data <- merged_data %>% filter(Hugo_Symbol %in% top_gene_names)

# Reshape data so that each row corresponds to a sample and each column corresponds to a gene
top_gene_data_wide <- top_gene_data %>%
  select(Sample_ID, Hugo_Symbol, Expression) %>%
  spread(key = Hugo_Symbol, value = Expression)

# Remove Sample_ID column for clustering
expression_matrix <- as.matrix(top_gene_data_wide[, -1])
rownames(expression_matrix) <- top_gene_data_wide$Sample_ID

# Determine optimal number of clusters using the elbow method
fviz_nbclust(expression_matrix, kmeans, method = "wss")

# Set the number of clusters (CURVE FLATENS AT POINT 3)
optimal_k <- 3

# Perform k-means clustering
set.seed(123) # For reproducibility
kmeans_result <- kmeans(expression_matrix, centers = optimal_k, nstart = 25)

# Add cluster assignments to the original data
top_gene_data_wide$Cluster <- kmeans_result$cluster

# Perform PCA
pca_result <- prcomp(expression_matrix, scale. = TRUE)

# Create a data frame with PCA results and cluster assignments
pca_data <- data.frame(Sample_ID = rownames(expression_matrix),
                       PC1 = pca_result$x[,1],
                       PC2 = pca_result$x[,2],
                       Cluster = as.factor(top_gene_data_wide$Cluster))

# Plot PCA
ggplot(pca_data, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(size = 3) +
  ggtitle("PCA of Top 100 Genes with K-means Clustering") +
  theme_minimal()


#-------------STEP 5: SURVIVAL ANALYSIS IN CLUSTERS----------------------------

# Randomly select one of the top 100 genes
set.seed(123)  # For reproducibility
selected_gene <- sample(top_genes$Hugo_Symbol, 1)
selected_gene

#----------Fit CoxPH Model for Each Cluster
# Load necessary library
library(survival)

# Extract necessary columns
survival_cols <- c("Sample_ID", "OS_Months", "OS_status")
gene_expr_col <- selected_gene

# Convert survival status to numeric (0 for censored, 1 for event)
merged_data$OS_status <- as.numeric(gsub(".*:.*", "\\1", merged_data$OS_status))

# Iterate through each cluster and fit CoxPH model
coxph_results <- list()

# Loop through each cluster
for (cluster in unique(top_gene_data_wide$Cluster)) {
  # Filter data for the current cluster
  cluster_sample_ids <- top_gene_data_wide %>% filter(Cluster == cluster) %>% pull(Sample_ID)
  cluster_data <- merged_data %>%
    filter(Sample_ID %in% cluster_sample_ids & Hugo_Symbol == selected_gene) %>%
    select(Sample_ID, OS_Months, OS_status, Expression)
  
  if (nrow(cluster_data) > 0) {
    # Fit CoxPH model
    coxph_model <- coxph(Surv(OS_Months, OS_status) ~ Expression, data = cluster_data)
    coxph_summary <- summary(coxph_model)
    
    # Store results
    coxph_results[[paste("Cluster", cluster)]] <- list(
      model = coxph_model,
      summary = coxph_summary
    )
  }
}

# Print CoxPH results for each cluster
for (cluster in names(coxph_results)) {
  cat("\nResults for", cluster, ":\n")
  print(coxph_results[[cluster]]$summary)
}

#------------- Kaplan-Meier Survival Analysis--------------
library(dplyr)
perform_kaplan_meier_analysis <- function(cluster_data, selected_gene) {
  # Debugging: Print the first few rows of cluster_data
  print(head(cluster_data))
  
  # Split the data into high and low expression groups based on the median expression
  median_expression <- median(cluster_data$Expression, na.rm = TRUE)
  cluster_data <- cluster_data %>%
    mutate(Expression_Group = ifelse(Expression >= median_expression, "High", "Low"))
  
  # Debugging: Print the first few rows after grouping
  print(head(cluster_data))
  
  # Fit the Kaplan-Meier survival curves
  km_fit <- survfit(Surv(OS_Months, OS_status) ~ Expression_Group, data = cluster_data)
  
  # Plot the Kaplan-Meier curves
  plot <- ggsurvplot(km_fit, data = cluster_data, pval = TRUE,
                     title = paste("Kaplan-Meier Survival Curves for", selected_gene, "in Cluster"),
                     xlab = "Time (months)", ylab = "Survival Probability",
                     legend.title = "Expression Group")
  
  # Perform the log-rank test
  logrank_test <- survdiff(Surv(OS_Months, OS_status) ~ Expression_Group, data = cluster_data)
  
  return(list(plot = plot, logrank_test = logrank_test))
}

km_results <- list()

for (cluster in unique(top_gene_data_wide$Cluster)) {
  cluster_sample_ids <- top_gene_data_wide %>%
    filter(Cluster == cluster) %>%
    pull(Sample_ID)
  
  cluster_data <- merged_data %>%
    filter(Sample_ID %in% cluster_sample_ids & Hugo_Symbol == selected_gene) %>%
    select(Sample_ID, OS_Months, OS_status, Expression)
  
  # Debugging: Print the first few rows of cluster_data for each cluster
  print(paste("Cluster", cluster))
  print(head(cluster_data))
  
  if (nrow(cluster_data) > 0) {
    km_result <- perform_kaplan_meier_analysis(cluster_data, selected_gene)
    km_results[[paste("Cluster", cluster)]] <- km_result
    
    # Print plot for each cluster
    print(km_result$plot)
  } else {
    print(paste("No data for Cluster", cluster))
  }
}



