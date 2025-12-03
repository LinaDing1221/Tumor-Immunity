if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "limma", "org.Hs.eg.db"))
install.packages(c("dplyr", "pheatmap", "tibble"))

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(limma)
library(pheatmap)
library(org.Hs.eg.db)
library(tibble)
query <- GDCquery(
  project = "TCGA-COAD",     
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "HTSeq - Counts"
)
GDCdownload(query)
coad_data <- GDCprepare(query)
sample_barcode <- colnames(coad_data)
sample_type <- substr(sample_barcode, 14, 15)
sample_group <- case_when(
  sample_type %in% paste0("0", 1:9) ~ "Tumor",  
  sample_type == "11" ~ "Normal",              
  TRUE ~ "Other"                                
)
coad_data <- coad_data[, sample_group != "Other"]
sample_group <- factor(sample_group[sample_group != "Other"], 
                       levels = c("Normal", "Tumor"))
ferroptosis_genes <- c(
  "ALOX12", "DUOX1", "DUOX2", "CYBB", "NOX4", "SLC7A11", 
  "IREB2", "KRAS", "RB1", "EMC1", "EMC2", "NFE2L2", 
  "CS", "LPCAT3", "FTH1", "HSPB1", "RPL8", "HSF1", 
  "GPX4", "TP53", "G6PD", "PGD", "SQSTM1", "CARS1", 
  "SLC40A1", "NQO1", "FDFT1", "ATP5MC3", "VDAC2", 
  "AKR1C2", "AKR1C3", "ACSF2", "HMOX1", "MT1G"
)
gene_info <- rowData(coad_data) %>% as.data.frame()
gene_info$symbol <- mapIds(
  org.Hs.eg.db, keys = gene_info$gene_id,
  column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"
)
ferro_ensembl <- gene_info$gene_id[gene_info$symbol %in% ferroptosis_genes]
coad_ferro <- coad_data[ferro_ensembl, ]
count_matrix <- assays(coad_ferro)$"HTSeq - Counts"
keep <- rowSums(count_matrix > 10) >= 0.5 * ncol(count_matrix)
coad_ferro <- coad_ferro[keep, ]
count_matrix <- count_matrix[keep, ]
design <- model.matrix(~sample_group) 
voom_obj <- voom(count_matrix, design, plot = FALSE)
fit <- lmFit(voom_obj, design) %>% eBayes()
de_results <- topTable(
  fit, coef = 2, adjust = "fdr", 
  number = Inf, lfc = 1, p.value = 0.05
)
ferro_de_ensembl <- rownames(de_results)
expr_zscore <- t(scale(t(voom_obj$E[ferro_de_ensembl, ]))) 
rownames(expr_zscore) <- gene_info$symbol[match(rownames(expr_zscore), gene_info$gene_id)]
sample_anno <- data.frame(Group = sample_group)
rownames(sample_anno) <- colnames(expr_zscore)
pheatmap(
  expr_zscore,
  annotation_col = sample_anno,      
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),   show_rownames = TRUE,               
  show_colnames = FALSE,             
  cluster_rows = TRUE,               
  cluster_cols = TRUE,                
  treeheight_row = 10,
  treeheight_col = 10,
  main = "Ferroptosis-related DEGs in TCGA COAD (Normal vs Tumor)"
)
