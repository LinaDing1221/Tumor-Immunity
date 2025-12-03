if (!require("glmnet")) install.packages("glmnet")
if (!require("readxl")) install.packages("readxl")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
library(glmnet)
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
setwd("D:\\biowolf\\bioR\\38.forest")
data <- read_excel("Lasso-Riskscore.xlsx", sheet = "lambda.min")
target_genes <- c("ALOX12", "FDFT1", "KRAS", "HSPB1", "NOX1")
x <- as.matrix(data[, target_genes])  
y <- factor(data$event, levels = c(0, 1))  
set.seed(123)  
cv.lasso <- cv.glmnet(
  x = x, 
  y = y, 
  family = "binomial", 
  alpha = 1,  
  nfolds = 10, 
  type.measure = "deviance"  
)
lambda_min <- cv.lasso$lambda.min  
lambda_1se <- cv.lasso$lambda.1se 
cat("lambda.min =", round(lambda_min, 4), "\n")
cat("lambda.1se =", round(lambda_1se, 4), "\n")
coef_df <- coef(cv.lasso, s = "all") %>%  # 提取所有lambda的系数
  as.matrix() %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(log_lambda = log(rownames(.))) %>%  # 添加log(lambda)列
  select(log_lambda, all_of(target_genes)) %>%  # 保留log(lambda)和5个基因
  pivot_longer(cols = -log_lambda, names_to = "gene", values_to = "coefficient")  # 
trajectory_plot <- ggplot(coef_df, aes(x = log_lambda, y = coefficient, color = gene, linetype = gene)) +
  geom_line(linewidth = 1.2) +  
  geom_vline(xintercept = log(lambda_min), color = "red", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = log(lambda_1se), color = "blue", linetype = "dashed", linewidth = 1) +
  annotate("text", x = log(lambda_min), y = max(coef_df$coefficient)*0.8, 
           label = paste0("lambda.min = ", round(lambda_min, 4)), 
           color = "red", hjust = -0.1, size = 4) +
  annotate("text", x = log(lambda_1se), y = max(coef_df$coefficient)*0.7, 
           label = paste0("lambda.1se = ", round(lambda_1se, 4)), 
           color = "blue", hjust = -0.1, size = 4) +
  labs(x = "log(λ)", y = "LASSO Coefficient", color = "Gene", linetype = "Gene") +
  ggtitle("LASSO Regression Variable Trajectory Analysis") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
ggsave("LASSO_variable_trajectory.pdf", trajectory_plot, width = 10, height = 6, dpi = 300)
coef_min <- coef(cv.lasso, s = lambda_min) %>% as.matrix() %>% as.data.frame()
coef_1se <- coef(cv.lasso, s = lambda_1se) %>% as.matrix() %>% as.data.frame()
coef_spectrum <- data.frame(
  gene = target_genes,
  coef_lambda_min = coef_min[target_genes, 1],
  coef_lambda_1se = coef_1se[target_genes, 1]
) %>% 
  pivot_longer(cols = -gene, names_to = "lambda_type", values_to = "coefficient") %>% 
  mutate(lambda_type = factor(lambda_type, 
                              levels = c("coef_lambda_min", "coef_lambda_1se"),
                              labels = c("lambda.min", "lambda.1se")))
spectrum_plot <- ggplot(coef_spectrum, aes(x = gene, y = coefficient, fill = lambda_type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = round(coefficient, 3)), 
            position = position_dodge(width = 0.8), 
            vjust = ifelse(coef_spectrum$coefficient > 0, -0.3, 1.3),
            size = 3.5) +
  labs(x = "Gene", y = "LASSO Coefficient", fill = "Lambda Type") +
  ggtitle("LASSO Regression Coefficient Spectrum") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave("LASSO_coefficient_spectrum.pdf", spectrum_plot, width = 9, height = 6, dpi = 300)
