library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
setwd("D:\\biowolf\\bioR\\KEGGcircos") 
pFilter <- 0.05 
rt <- read.table("input.txt", sep="\t", header=T, check.names=F)
genes <- as.vector(rt[,1])  
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA) 
entrezIDs <- as.character(entrezIDs)
rt <- cbind(rt, entrezID=entrezIDs)
rt <- rt[is.na(rt[,"entrezID"])==F,]  
gene <- rt$entrezID  
go_bp <- enrichGO(
  gene = gene,
  OrgDb = org.Hs.eg.db,  
  ont = "BP",  
  keyType = "ENTREZID",  
  pvalueCutoff = pFilter,  
  qvalueCutoff = 1,  
  readable = TRUE  
)
go_result <- as.data.frame(go_bp)
if (nrow(go_result) == 0) stop
go_result$GeneRatio <- sapply(
  strsplit(go_result$GeneRatio, "/"),
  function(x) as.integer(x[1]) / as.integer(x[2])
)
go_result <- go_result[order(-go_result$GeneRatio), ]
go_result$Term <- factor(go_result$Term, levels = go_result$Term)
pdf(file="GO_BP_bubble.pdf", width=8, height=6)  
ggplot(go_result, aes(x = GeneRatio, y = Term)) +
  geom_point(aes(size = Count, color = pvalue), alpha = 0.9) +
  scale_color_gradient(
    low = "#0088FF", 
    high = "#FF5555", 
    name = "P value",
    breaks = c(0.01, 0.02, 0.03),
    labels = c("0.01", "0.02", "0.03")
  ) +
  scale_size_continuous(
    name = "Counts",
    breaks = c(1, 2),  
    range = c(4, 8),  
    labels = c("1.00", "2.00")
  ) +
  labs(x = "GeneRatio", y = "") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 10, color = "black"),
    legend.text = element_text(size = 9, color = "black"),
    legend.box = "vertical",  
    plot.margin = margin(10, 20, 10, 10, "pt")  
  )
dev.off()
