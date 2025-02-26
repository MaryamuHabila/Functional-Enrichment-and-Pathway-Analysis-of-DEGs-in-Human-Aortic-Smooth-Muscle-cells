# Functional-Enrichment-and-Pathway-Analysis-of-DEGs-in-Human-Aortic-Smooth-Muscle-cells
This project aims to analyze differentially expressed genes (DEGs) in Human Aortic Smooth Muscle Cells (HAoSMCs) using functional enrichment and pathway analysis techniques. The primary objective is to link upregulated genes to biological functions, molecular mechanisms, and metabolic/signaling pathways. The analysis integrates Gene Ontology (GO) and pathway enrichment tools to provide a comprehensive understanding of gene regulation in HAoSMCs.

Objectives:

Data Preparation: Import HAoSMCs differential expression results from an Excel file. Ensure log2 fold change and p-values are correctly formatted as numeric values. Perform data cleaning to remove non-numeric characters and missing values. Compute -log10(p-value) for enhanced visualization of statistical significance.

Upregulated Gene Selection: Identify significantly upregulated genes with a log2 fold change > 1.5. Extract gene identifiers for functional enrichment analysis.

Gene Ontology (GO) Enrichment Analysis: Perform GO analysis using the clusterProfiler package in R. Categorize results into:

Biological Processes (BP): Investigates roles in physiological functions.

Molecular Functions (MF): Examines molecular activities and interactions.

Cellular Components (CC): Identifies subcellular locations of gene products.

Visualize the top 50 most significant GO terms using bar plots.

Pathway Analysis: Investigate metabolic and signaling pathways using KEGG. Map DEGs to curated biological pathways using ReactomePA. Explore disease-related pathways and their potential implications.

Data Visualization: Generate bar plots to highlight key GO terms in BP, MF, and CC categories. Rank GO terms by -log10(p-value) to prioritize the most significant biological insights.
Save high-resolution visualizations for reporting and presentation.

Tools and Libraries Used: R packages: readxl, ggplot2, dplyr, ggrepel, writexl, ggraph, clusterProfiler, org.Hs.eg.db, ReactomePA

Data Source: HAoSMCs differential expression results processed via DESeq2

Significance of the Study: This project provides a systematic approach to understanding the biological significance of upregulated genes in HAoSMCs. The insights gained from this analysis can help in: Identifying key regulatory mechanisms in smooth muscle cell function. Understanding potential molecular targets in cardiovascular diseases. Supporting future experimental validation of key pathways. By integrating GO enrichment and pathway analysis, this study bridges the gap between gene expression changes and functional relevance in cellular and disease contexts.


