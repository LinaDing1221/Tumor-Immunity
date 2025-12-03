setwd("D:\\biowolf\\bioR\\38.forest")
rt <- read.table("input.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
rt <- rt[order(-rt$HR), ]
gene <- rownames(rt)  
n <- nrow(rt)        
hr <- sprintf("%.3f", rt$HR)
hrLow <- sprintf("%.3f", rt$HR.95L)
hrHigh <- sprintf("%.3f", rt$HR.95H)
Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")") 
pVal <- rt$pvalue
pVal_label <- ifelse(pVal < 0.001, "***",
                     ifelse(pVal < 0.01, "**",
                            ifelse(pVal < 0.05, "*", sprintf("%.3f", pVal))))
pVal_text <- ifelse(pVal < 0.001, "<0.001", sprintf("%.3f", pVal))
pVal_final <- paste0(pVal_text, " ", pVal_label)  
pdf(file = "forest.pdf", width = 8, height = 10)  
nRow <- n + 1  
ylim <- c(1, nRow)
layout(matrix(c(1, 2), nc = 2), width = c(4, 4))
par(mar = c(4, 2, 2, 1.5)) 
xlim_text <- c(0, 5)
plot(1, xlim = xlim_text, ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
text.cex <- ifelse(n > 25, 0.65, 0.75)
text(0.1, n:1, gene, adj = 0, cex = text.cex, font = 1)
text(2.2, n:1, pVal_final, adj = 1, cex = text.cex)
text(2.2, n + 1, "P value", cex = text.cex + 0.1, font = 2, adj = 1)  
text(4.8, n:1, Hazard.ratio, adj = 1, cex = text.cex)
text(4.8, n + 1, "Hazard Ratio (95%CI)", cex = text.cex + 0.1, font = 2, adj = 1)  
par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))  
x_min <- min(as.numeric(hrLow)) * 0.9
x_max <- max(as.numeric(hrHigh)) * 1.1
xlim_forest <- c(x_min, x_max)
plot(1, xlim = xlim_forest, ylim = ylim, type = "n", axes = FALSE, ylab = "", 
     xaxs = "i", xlab = "Hazard Ratio")
arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, 
       angle = 90, code = 3, length = 0.025, col = "darkblue", lwd = 2)
boxcolor <- ifelse(as.numeric(hr) > 1, "red", 
                   ifelse(as.numeric(hr) < 1, "blue", "black"))
points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.2)
abline(v = 1, col = "black", lty = 2, lwd = 2)
axis(1, at = seq(round(x_min, 1), round(x_max, 1), by = 0.2), cex.axis = text.cex)
title(main = "Univariate Cox Regression Analysis of Genes", 
      line = 0.5, cex.main = 1.2, font.main = 2)
dev.off()
