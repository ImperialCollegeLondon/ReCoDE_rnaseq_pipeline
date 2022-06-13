[1mdiff --git a/notebooks/downstream_analysis.Rmd b/notebooks/downstream_analysis.Rmd[m
[1mindex f594f67..0608acb 100644[m
[1m--- a/notebooks/downstream_analysis.Rmd[m
[1m+++ b/notebooks/downstream_analysis.Rmd[m
[36m@@ -18,7 +18,7 @@[m [mlibrary(edgeR)  # https://bioconductor.org/packages/release/bioc/html/edgeR.html[m
 library(SummarizedExperiment)  # https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html[m
 [m
 # package for plotting PCA[m
[31m-library(ggbiplot)  # https://github.com/vqv/ggbiplot[m
[32m+[m[32mlibrary(ggplot2)  #[m[41m [m
 [m
 # package for making volcano plots[m
 library(EnhancedVolcano)  # https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html[m
[36m@@ -168,7 +168,7 @@[m [msample groups.[m
 ```{r}[m
 [m
 # select genes to keep[m
[31m-fbe_keep <- filterByExpr(se, group = se$condition)[m
[32m+[m[32mfbe_keep <- filterByExpr(se, group = se$condition, min.count = 500)[m[41m[m
 [m
 # remove genes that are not in the list[m
 filtered_se <- se[fbe_keep, ][m
[36m@@ -246,7 +246,7 @@[m [mstandardisation is recommended to keep different genes on comparable scales.[m
 # prcomp calculates PCA using singular value decomposition[m
 cpm_pca <- prcomp(t(assay(filtered_se, 2)), center = TRUE, scale. = TRUE)[m
 [m
[31m-print(cpm_pca)[m
[32m+[m[32mprint(cpm_pca$x)[m[41m[m
 [m
 ```[m
 [m
[36m@@ -263,10 +263,11 @@[m [mwithout knowledge of labels of interest.[m
 [m
 ```{r}[m
 [m
[31m-# use the ggbiplot package to plot the PCA[m
[31m-ggbiplot(cpm_pca, groups = filtered_se$condition, ellipse = TRUE) +[m
[31m-  scale_color_discrete(name = '') +[m
[31m-  theme(legend.position = 'right')[m
[32m+[m[32mpca_df <- data.frame(condition = filtered_se$condition, cpm_pca$x)[m[41m[m
[32m+[m[41m[m
[32m+[m[32m# use the ggplot package to plot the PCA[m[41m[m
[32m+[m[32mggplot(pca_df, aes(PC1, PC2, col = condition)) +[m[41m[m
[32m+[m[32m    geom_point()[m[41m[m
 [m
 ```[m
 [m
