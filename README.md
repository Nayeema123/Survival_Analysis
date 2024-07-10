# TCGA-BRCA Dataset Analysis: Gene Expression and Survival in Breast Cancer

This project explores the relationship between gene expression profiles and survival outcomes in breast cancer using data from The Cancer Genome Atlas (TCGA) Breast Invasive Carcinoma (BRCA) dataset. The analysis integrates genomic data with clinical information to identify genes that significantly impact patient survival and to understand how these genes manifest in distinct patient clusters.

# Objectives
# 1. Data Acquisition: 
   > Downloaded gene expression and clinical data from cBioPortal to ensure comprehensive coverage of patient genomic profiles and clinical outcomes.

# 2. Data Preparation:
   > Preprocessing: Filtered gene expression data to remove noise and irrelevant features, ensuring high-quality data for subsequent analysis.

   > Survival Data Extraction: Extracted survival data, including time-to-event and event status, from clinical records to enable survival analysis.

# 3. Feature Selection:
   > Top Gene Identification: Employed statistical methods such as Cox proportional hazards models and the log-rank test to identify the top 100 genes significantly associated with 
     patient survival. These genes serve as the basis for further analysis.

# 4. Clustering analysis:
   > K-means Clustering: Applied k-means clustering on the selected genes to categorize patients into clusters based on their gene expression profiles.
   
   > Cluster Validation: Determined the optimal number of clusters using techniques like the elbow method, ensuring robust cluster assignments.
   
# 5. Survival analysis:
   > Within-cluster Analysis: Conducted Cox proportional hazards modeling within each cluster to assess how gene expression influences survival outcomes.

   > Kaplan-Meier Curves: Plotted Kaplan-Meier survival curves for high and low expression groups of a representative gene in each cluster.

   > Statistical Validation: Performed log-rank tests to statistically compare survival distributions between these groups, evaluating the significance of gene expression levels on 
     patient outcomes.
   
# Findings and Insights
> Identified specific genes that significantly affect survival patterns across different patient clusters.

> Characterized distinct survival profiles associated with high and low expression levels of key genes within each cluster.

> Provided insights into potential biomarkers and therapeutic targets based on gene expression patterns and survival associations.
