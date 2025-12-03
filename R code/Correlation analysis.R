if (!require("corrplot")) install.packages("corrplot")
if (!require("igraph")) install.packages("igraph")
if (!require("reshape2")) install.packages("reshape2")
library(corrplot)
library(igraph)
library(reshape2)

setwd("D:\\biowolf\\bioR\\23.corrplot")
iron_genes <- c("ALOX12", "FDFT1", "KRAS", "HSPB1", "NOX1")  
sting_genes <- c("STING1", "TBK1", "IRF3", "MAP3K6", "QSOX1", "SLC35A4")  
rt <- read.table("input.txt", sep="\t", header=T, row.names=1, check.names=F)
rt_sub <- rt[, c(iron_genes, sting_genes)]  
x_iron <- rt_sub[, iron_genes] 
x_sting <- rt_sub[, sting_genes]  
cor_mat <- cor(x_iron, x_sting)
calc_cor_p <- function(x, y) {
  p_mat <- matrix(NA, nrow = ncol(x), ncol = ncol(y))
  rownames(p_mat) <- colnames(x)
  colnames(p_mat) <- colnames(y)
  for (i in 1:ncol(x)) {
    for (j in 1:ncol(y)) {
      p_mat[i, j] <- cor.test(x[, i], y[, j])$p.value
    }
  }
  return(p_mat)
}
p_mat <- calc_cor_p(x_iron, x_sting)
pdf(file="STING_iron_corr.pdf", width=8, height=6)
sig_label <- matrix("", nrow=nrow(cor_mat), ncol=ncol(cor_mat))
sig_label[p_mat < 0.001] <- "***"
sig_label[p_mat >= 0.001 & p_mat < 0.01] <- "**"
sig_label[p_mat >= 0.01 & p_mat < 0.05] <- "*"
corrplot(
  cor_mat,
  method = "color",         
  col = colorRampPalette(c("blue", "white", "red"))(50),  
  tl.pos = "lt",             
  tl.col = "black",         
  tl.cex = 1.0,              
  cl.pos = "r",              
  cl.cex = 0.9,              
  addgrid.col = NA,          
  is.corr = TRUE,            
  mar = c(2, 4, 4, 2)        
)
text(
  x = col(cor_mat), 
  y = nrow(cor_mat) - row(cor_mat) + 1,  
  labels = sig_label, 
  cex = 1.2, 
  col = "black"
)
legend(
  "bottomright", 
  legend = c("*** p<0.001", "** p<0.01", "* p<0.05"),
  bty = "n", 
  cex = 0.9
)
dev.off()
outFile_net <- "STING_iron_corNetwork.pdf"
cutoff <- 0.4  
full_cor <- cor(rt_sub)
full_cor[upper.tri(full_cor)] <- NA  
df_melt <- melt(as.data.frame(full_cor), id.vars = rownames(full_cor), variable.name = "gene2", value.name = "cor")
df_melt <- df_melt[!is.na(df_melt$cor) & df_melt$Var1 != df_melt$gene2 & abs(df_melt$cor) > cutoff, ]
g <- graph.data.frame(df_melt[, c("Var1", "gene2")], directed = FALSE)
E(g)$color <- ifelse(df_melt$cor > 0, rgb(254,67,101, maxColorValue=255, alpha=abs(df_melt$cor)*255), 
                     rgb(0,0,255, maxColorValue=255, alpha=abs(df_melt$cor)*255))
V(g)$size <- 10
V(g)$label.cex <- 1.1
V(g)$color <- "white"
V(g)$frame.color <- "black"
E(g)$width <- abs(df_melt$cor)*5 
pdf(outFile_net, width=9, height=8)
layout(matrix(c(1,1,0,2,0), byrow=T, nc=3), height=c(8,1), width=c(4,1,4))
par(mar=c(2,2,2,2))
plot(g, layout=layout_nicely, edge.arrow.size=0, vertex.label.color="black", edge.curved=0.3)
par(mar=c(2,1,1,1), xpd=TRUE)
color_legend <- c(rgb(254,67,101, seq(255,0,by=-2.55), maxColorValue=255),
                  rgb(0,0,255, seq(0,255,by=2.55), maxColorValue=255))
barplot(rep(1, length(color_legend)), border=NA, space=0, axes=F, col=color_legend)
axis(3, at=seq(1, length(color_legend), length=5), labels=c("1.0", "0.5", "0", "-0.5", "-1.0"), tick=F)
dev.off()

