#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("ReactomePA")
#####REACTOME 

############ PATHWAY ANALYSIS ###############
#This script performs GO enrichment analysis on a list of upregulated genes, sorts the most significant terms, and visualizes the top 50 using a bar plot.
#Identified biological functions, Molecular Functions and cellular components
#Explored metabolic and signaling pathways (KEGG)
#Mapped genes to curated biological pathways (Reactome)
#In summary, this analysis helps understand the functional roles of differentially expressed genes in a biological system by linking them to relevant biological processes (BP), molecular functions (MF), cellular components (CC), metabolic pathways, and disease-related pathways.


# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(writexl)
library(ggraph)
library(clusterProfiler) #Provides functions for functional enrichment analysis, including KEGG.
library(org.Hs.eg.db) #A gene annotation database for humans, mapping gene symbols to identifiers.
library(ReactomePA)


#### IMPORTING HAoSMCs Excel Data ####
#Specify the path and read the Excel file containing HAoSMCs (Human Aortic Smooth Muscle Cells) DESeq2 results into R. It specifies the file path and uses the read_excel() function to load the data into a variable (DESeq_HAoSMCs_subset4). 
file_path_subset4 <- "/Users/maryamu/Desktop/SMCs_DESeq2_results_subset(n=3).xlsx"
# Read the Excel file (default is the first sheet)
DESeq_HAoSMCs_subset4 <- read_excel(file_path_subset4)


#Convert log2FoldChange to Numeric
#Converts the log2FoldChange column from character/factor to numeric format.If there are non-numeric values, they will be converted to NA
DESeq_HAoSMCs_subset4$log2FoldChange <- as.numeric(DESeq_HAoSMCs_subset4$log2FoldChange)


######. Ensure log2FoldChange and pvalue are numeric. #####
#Check for Non-Numeric Values in log2FoldChange.Displays the unique values in the column.Helps identify if there are any unexpected values (e.g., text, symbols).
unique(DESeq_HAoSMCs_subset4$log2FoldChange)


#Identify Problematic Rows (Non-Numeric Values) in log2FoldChange.Finds rows where conversion to numeric failed (NA values).Useful for debugging issues with data formatting.
DESeq_HAoSMCs_subset4[is.na(as.numeric(DESeq_HAoSMCs_subset4$log2FoldChange)), ]


#Remove Non-Numeric Characters in log2FoldChange.Uses gsub() to remove unwanted characters (anything that isn’t a number, decimal, or negative sign). This is important if there are spaces, special symbols, or text values that prevent conversion to numeric.
DESeq_HAoSMCs_subset4$log2FoldChange <- gsub("[^0-9.-]", "", DESeq_HAoSMCs_subset4$log2FoldChange)



#Convert Cleaned log2FoldChange to Numeric
#Ensures that log2FoldChange is properly formatted as numeric after cleaning.
DESeq_HAoSMCs_subset4$log2FoldChange <- as.numeric(DESeq_HAoSMCs_subset4$log2FoldChange)



#Remove Rows with NA in log2FoldChange
#Excludes any rows where log2FoldChange is still NA (invalid values)
DESeq_HAoSMCs_subset4 <- DESeq_HAoSMCs_subset4[!is.na(DESeq_HAoSMCs_subset4$log2FoldChange), ]


# Ensure log2FoldChange are numeric
DESeq_HAoSMCs_subset4$log2FoldChange <- as.numeric(DESeq_HAoSMCs_subset4$log2FoldChange)


#####Repeat the Process for p-value.The exact same steps are applied to the p-value column to ensure it is clean and numeric.
DESeq_HAoSMCs_subset4$pvalue <- as.numeric(DESeq_HAoSMCs_subset4$pvalue)
unique(DESeq_HAoSMCs_subset4$pvalue)
#Identify Problematic Rows: Locate rows where the conversion to numeric fails.
DESeq_HAoSMCs_subset4[is.na(as.numeric(DESeq_HAoSMCs_subset4$pvalue)), ]
#Remove Non-Numeric Characters
DESeq_HAoSMCs_subset4$pvalue <- gsub("[^0-9.-]", "", DESeq_HAoSMCs_subset4$pvalue)
DESeq_HAoSMCs_subset4$pvalue <- as.numeric(DESeq_HAoSMCs_subset4$pvalue)
#Exclude rows with NA values:
DESeq_HAoSMCs_subset4 <- DESeq_HAoSMCs_subset4[!is.na(DESeq_HAoSMCs_subset4$pvalue), ]
# Ensure pvalue are numeric
DESeq_HAoSMCs_subset4$pvalue <- as.numeric(DESeq_HAoSMCs_subset4$pvalue)

###### Calculate -log10(pvalue). #######
#Computes the negative log10 transformation of the p-value, which is commonly used in volcano plots and significance testing.This transformation helps visualize small p-values more effectively
DESeq_HAoSMCs_subset4$neg_log10_pvalue <- -log10(DESeq_HAoSMCs_subset4$pvalue)
#Obtain the Length of the Dataset
length(DESeq_HAoSMCs_subset4) #Returns the number of columns in the dataset,
nrow(DESeq_HAoSMCs_subset4 #Returns the number of rows


##### SUMMARY ######
#Cleans log2FoldChange and p-value by removing non-numeric characters.
#Ensures these columns are properly formatted as numeric values.
#Removes invalid rows where conversion failed.
#Calculates -log10(p-value) for downstream statistical analysis.
#This ensures that the dataset is ready for differential expression analysis and visualization.



########### PATHWAY ANALYSES ###############

#### UREGULATED ##########

# Step 1: Extract upregulated genes (Log2FC > 1.5)
#This step filters the dataset to keep only upregulated genes, defined as those with a log2 fold change (Log2FC) greater than 1.5.filter(log2FoldChange > 1.5) → Selects rows where log2FoldChange is greater than 1.5, meaning these genes are significantly upregulated. %>% (Pipe Operator) → Passes DESeq_HAoSMCs_subset4 into filter(), making the code more readable.upregulated_genes_filtered → Stores the filtered dataset for further analysis.This step helps in identifying genes that show a strong increase in expression in the experimental condition compared to the control.
upregulated_genes_filtered <- DESeq_HAoSMCs_subset4 %>%
  filter(log2FoldChange > 1.5)


###Gene Ontology (GO) Enrichment Analysis and Visualization
#This script performs GO enrichment analysis on a list of upregulated genes and visualizes the top 50 enriched GO terms using a bar plot.
#Performs GO enrichment analysis using the enrichGO function from clusterProfiler
#The analysis focuses on Biological Processes (BP), Molecular Function (MF) and Cellular Component (CC)
#No filtering on p-value/q-value so that all GO terms are included initially.
#Extract gene IDs for GO analysis
gene_list_filtered <- upregulated_genes_filtered$gene_id  # Extracts the gene_id column from the filtered list of upregulated genes.Ensures the column name matches the required gene identifier type (e.g., "SYMBOL" or "ENTREZID").

# Step 3: Run GO Enrichment Analysis without filtering on p-value or q-value. Select Top 50 GO Terms for visualisation.Sorts the GO results by p-value to prioritize the most statistically significant terms.
go_results <- enrichGO(
  gene          = gene_list_filtered,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  # Adjust to your gene ID type (e.g., "ENTREZID" or "SYMBOL")
  ont           = "BP",     # Options: "BP", "CC", "MF", or "ALL" for all categories
  pAdjustMethod = "BH",      # Benjamini-Hochberg correction
  pvalueCutoff  = 1,         # Include all results regardless of p-value
  qvalueCutoff  = 1          # Include all results regardless of q-value
)

# Step 4: Select the top 50 terms based on p-value. Sorts the GO results by p-value to prioritize the most statistically significant terms.
top_go_results <- go_results[order(go_results@result$pvalue), ]  # Sort by p-value
top_50_go_results <- head(top_go_results, 50)  # Select top 50 terms

# Step 5: Visualize the top 50 GO terms. Creates a data frame for visualization.Converts p-values to -log10(p-value) to make significant terms stand out. 
barplot_data <- data.frame(
  Description = top_50_go_results$Description,
  neg_log10_pvalue = -log10(top_50_go_results$pvalue)
)

# Ensure data is sorted by significance.Sorts terms in descending order of significance.
barplot_data <- barplot_data[order(barplot_data$neg_log10_pvalue, decreasing = TRUE), ]

# Bar plot of the top 50 GO terms
# Ensure the data frame is sorted by -log10(p-value) in descending order
barplot_data <- barplot_data %>%
  arrange(desc(neg_log10_pvalue))

# Create a bar plot with the longest bars at the top.Generate a Bar Plot
#Uses ggplot2 to create a horizontal bar plot of the top 50 GO terms.GO terms are sorted by significance (-log10 p-value).The most significant terms appear at the top.
ggplot(barplot_data, aes(x = reorder(Description, neg_log10_pvalue), y = neg_log10_pvalue)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip coordinates for horizontal bars
  labs(
    title = "Top 50 GO Terms (Ordered by -log10(p-value))",
    x = "GO Terms",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 13))

# Save the Bar Plot.Saves the plot as a high-resolution PNG file (300 dpi) for better visualization.
ggsave(
  filename = "~/Documents/GO_Terms_Barplot.png", 
  plot = last_plot(),  # Saves the most recently created ggplot object
  width = 10, height = 6, dpi = 300
)



###Cellular component
#GO Enrichment Analysis for Cellular Components (CC) and Visualization
#This script performs Gene Ontology (GO) enrichment analysis for the Cellular Component (CC) ontology, identifies the top 50 most significant GO terms, and visualizes them using a bar plot.
#Run GO Enrichment Analysis without filtering on p-value or q-value
#Runs GO enrichment analysis using the enrichGO function from clusterProfiler.
#Focuses on Cellular Component (CC) ontology, meaning it finds enriched terms related to cellular structures (e.g., nucleus, membrane, ribosome).
#Does not filter by p-value or q-value, allowing all results to be included
go_results <- enrichGO(
  gene          = gene_list_filtered, #List of genes for analysis
  OrgDb         = org.Hs.eg.db, #Database for human genes
  keyType       = "SYMBOL",  # Adjust to your gene ID type (e.g., "ENTREZID" or "SYMBOL")
  ont           = "CC",     # Options: "BP", "CC", "MF", or "ALL" for all categories
  pAdjustMethod = "BH",      # Adjusts Benjamini-Hochberg correction
  pvalueCutoff  = 1,         # Include all results regardless of p-value
  qvalueCutoff  = 1          # Include all results regardless of q-value
)

#Select the top 50 terms based on p-value or other criteria. #Sorts GO terms by p-value, prioritizing the most statistically significant ones.
top_go_results <- go_results[order(go_results@result$pvalue), ]  # Sort by p-value
top_50_go_results <- head(top_go_results, 50)  # Select top 50 terms

# Extracts the top 50 GO terms for visualization.Creates a data frame with GO term descriptions and their -log10(p-value) values. # Ensure data is sorted by significance. Sorting by -log10(p-value) ensures the most significant terms appear at the top of the plot.
barplot_data <- data.frame(
  Description = top_50_go_results$Description,
  neg_log10_pvalue = -log10(top_50_go_results$pvalue)
)
barplot_data <- barplot_data[order(barplot_data$neg_log10_pvalue, decreasing = TRUE), ] 

# Bar plot of the top 50 GO terms
# Ensure the data frame is sorted by -log10(p-value) in descending order
barplot_data <- barplot_data %>%
  arrange(desc(neg_log10_pvalue))

#Creates a horizontal bar plot using ggplot2.
#GO terms are ordered by -log10(p-value) so that the most significant ones are at the top.
#Uses a minimal theme for better readability.
# Create a bar plot with the longest bars at the top
ggplot(barplot_data, aes(x = reorder(Description, neg_log10_pvalue), y = neg_log10_pvalue)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip coordinates for horizontal bars
  labs(
    title = "Top 50 GO Terms (Ordered by -log10(p-value))",
    x = "GO Terms",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8))

# Save the Bar Plot
#Saves the plot as a high-resolution PNG file (300 dpi) for better visualization
ggsave(
  filename = "~/Documents/GO_Terms_Barplot.png", 
  plot = last_plot(),  # Saves the most recently created ggplot object
  width = 10, height = 6, dpi = 300
)



###Molecular function
# Step 3: Run GO Enrichment Analysis without filtering on p-value or q-value
go_results <- enrichGO(
  gene          = gene_list_filtered,
  OrgDb         = org.Hs.eg.db,
  keyType       = "SYMBOL",  # Adjust to your gene ID type (e.g., "ENTREZID" or "SYMBOL")
  ont           = "MF",     # Options: "BP", "CC", "MF", or "ALL" for all categories
  pAdjustMethod = "BH",      # Benjamini-Hochberg correction
  pvalueCutoff  = 1,         # Include all results regardless of p-value
  qvalueCutoff  = 1          # Include all results regardless of q-value
)

# Step 4: Select the top 50 terms based on p-value or other criteria
top_go_results <- go_results[order(go_results@result$pvalue), ]  # Sort by p-value
top_50_go_results <- head(top_go_results, 50)  # Select top 50 terms

# Step 5: Visualize the top 50 GO terms
# Create a bar plot
barplot_data <- data.frame(
  Description = top_50_go_results$Description,
  neg_log10_pvalue = -log10(top_50_go_results$pvalue)
)

# Ensure data is sorted by significance
barplot_data <- barplot_data[order(barplot_data$neg_log10_pvalue, decreasing = TRUE), ]

# Bar plot of the top 50 GO terms

# Ensure the data frame is sorted by -log10(p-value) in descending order
barplot_data <- barplot_data %>%
  arrange(desc(neg_log10_pvalue))

# Create a bar plot with the longest bars at the top
ggplot(barplot_data, aes(x = reorder(Description, neg_log10_pvalue), y = neg_log10_pvalue)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip coordinates for horizontal bars
  labs(
    title = "Top 50 GO Terms (Ordered by -log10(p-value))",
    x = "GO Terms",
    y = "-log10(p-value)"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(size = 8))





########. KEGG enrichment ######
#KEGG (Kyoto Encyclopedia of Genes and Genomes) enrichment analysis identifies biological pathways significantly associated with a set of genes. It helps understand how genes interact in functional pathways, such as metabolism, disease mechanisms, and signaling networks.
# Convert gene symbols to Entrez IDs
#Many enrichment tools (like KEGG) use Entrez Gene IDs instead of gene symbols.
#bitr() converts gene symbols to Entrez IDs using the org.Hs.eg.db database.
#The converted Entrez IDs are stored in gene_list_entrez.
gene_list_entrez <- bitr(gene_list_filtered, 
                         fromType = "SYMBOL", 
                         toType = "ENTREZID", 
                         OrgDb = org.Hs.eg.db)
# Use the converted Entrez IDs
gene_list_entrez <- gene_list_entrez$ENTREZID
#Perform KEGG Enrichment Analysis
#Identifies significantly enriched KEGG pathways associated with the input genes.
#The organism = "hsa" option specifies that the analysis is for human genes.
#No p-value filtering (pvalueCutoff = 1) ensures that all results are included.
kegg_results <- enrichKEGG(gene = gene_list_entrez, 
                           organism = "hsa", # "hsa" stands for Homo sapiens (human)
                           pvalueCutoff = 1 #Includes all results regardless of p-value)

#  Visualize KEGG Enrichment Results with a Dot Plot
# Create dot plot with larger, bold KEGG terms
#Creates a dot plot of the top 30 enriched KEGG pathways.The size of dots represents the number of genes involved in each pathway.The color of dots represents the statistical significance (p-value).

p <- dotplot(kegg_results, showCategory = 30) +
  ggtitle("KEGG Pathway Enrichment") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(face = "bold", size = 27),  # Larger & bold KEGG terms
    axis.title.y = element_text(size = 18, face = "bold"), # Bigger y-axis title
    axis.text.x = element_text(size = 14),                 # Increase x-axis text size
    axis.title.x = element_text(size = 16, face = "bold"), # Bigger x-axis title
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)  # Bigger, centered title
  )

# Save as high-resolution 300dpi image
ggsave("KEGG_Enrichment.png", plot = p, dpi = 300, width = 12, height = 7, units = "in")


######## Reactome Pathway Enrichment Analysis ######
#Reactome pathway enrichment analysis identifies biological pathways associated with a given set of genes. It helps understand how genes function together in processes like cell signaling, metabolism, and disease mechanisms.
#Step 1: Perform Reactome Pathway Enrichment Analysis
#enrichPathway() is a function from the ReactomePA package that identifies pathways enriched in the input gene list.
#gene = gene_list_entrez: Uses Entrez Gene IDs as input for enrichment.
#organism = "human": Specifies that the analysis is for Homo sapiens.
#pvalueCutoff = 1: Includes all pathways without filtering based on significance.
#Creates a dot plot of the top 15 enriched Reactome pathways. Dot size represents the number of genes involved in each pathway.Dot color represents the p-value (significance)

reactome_results <- enrichPathway(gene = gene_list_entrez, 
                                  organism = "human", 
                                  pvalueCutoff = 1 ## Include all results regardless of p-value)
                                  
#Visualize the Top 15 Enriched Pathways Using a Dot Plot. 
dotplot(reactome_results, showCategory = 15)

#Save as a High-Resolution (300 dpi) Image
ggsave(
  "/Users/maryamu/Documents/reactome_smcsupregpvaldotplot1.png", 
  plot = dotplot(reactome_results, showCategory = 20), 
  width = 12, height = 10, dpi = 300
)




