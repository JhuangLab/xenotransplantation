.libPaths(new = c(.libPaths(), "/cluster/home/danyang_jh/sbin/R/R-4.3.0", "/cluster/home/danyang_jh/sbin/R/R-4.2.1"))
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "QFeatures", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "MSnbase", "msqrob2", 
          "enrichplot", "org.Hs.eg.db", "parallel", "jhuanglabRNAseq", "jhuanglabGO", "limma", 
          "readxl", "writexl", "ComplexHeatmap", "circlize", "clusterProfiler", "MsCoreUtils", 
          "missForest", "paletteer", "factoextra", "FactoMineR", "openxlsx","gdata","ggrepel")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "liver"
dataset <- "wangrongrong"
species <- "Xeno_pig_monkey"
workdir <- glue::glue("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA")
workdir %>% fs::dir_create() %>% setwd()


output_dir <- "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA"

file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_mmul.csv"
expr_matrix <- read_csv(file_path)

picture_list <- list()

results <- expr_matrix
comp <- "LT_vs_LC_mmul"
results <- results %>% as.data.frame() %>% 
  mutate(trend = case_when((PValue < 0.01) & (logFC > 0.5) ~ "up", 
                           (PValue < 0.01) & (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
results <- results[!duplicated(results$gene_name),]
top_diff_genes <- results[order(results$PValue), ][1:10, ] %>% .[.$trend != "Not-sig", "gene_name"]

label_frame <- data_frame(gene_name = c(unique(top_diff_genes)),
                          tag = c(unique(top_diff_genes)))
results <- left_join(results, label_frame, by = "gene_name")
## volcanic
vol1 <- ggplot(results, aes(x = logFC, y = -log10(PValue), color = trend, fill = trend)) +
  ylim(0, max(-log10(results$PValue))+0.5)+
  xlim(-35,35)+
  geom_point(cex = 1) + 
  ggrepel::geom_label_repel(aes(label = tag), 
                            size = 2, 
                            box.padding = 0.5, 
                            max.overlaps = 40, 
                            show.legend = FALSE, 
                            color = "white",  # 设置文本颜色为白色
                            label.size = 0,
                            segment.color = "black")+  # 去掉标签边框
  #ggrepel::geom_text_repel(aes(label = tag), size = 2, box.padding = 0.5, max.overlaps = 40,show.legend = FALSE) +
  scale_color_manual(values = alpha(c("Not-sig" = "gray50", "up" = "red", "down" = "blue"), 0.75)) +
  scale_fill_manual(values = c("Not-sig" = "gray50", "up" = "red", "down" = "blue")) +  # 设置填充颜色
  
  my_theme1 + labs(title = comp)


ggsave(glue(output_dir,"/single_volcano_filtered_pval_mmul_V2.pdf"), vol1, width = 5, height = 3)



results <- expr_matrix
comp <- "LT_vs_LC_mmul"
results <- results %>% as.data.frame() %>% 
  mutate(trend = case_when((fdr < 0.05) & (logFC > 0.5) ~ "up", 
                           (fdr < 0.05) & (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
results <- results[!duplicated(results$gene_name),]
top_diff_genes <- results[order(results$fdr), ][1:10, ] %>% .[.$trend != "Not-sig", "gene_name"]

label_frame <- data_frame(gene_name = c(unique(top_diff_genes)),
                          tag = c(unique(top_diff_genes)))
results <- left_join(results, label_frame, by = "gene_name")
## volcanic
vol1 <- ggplot(results, aes(x = logFC, y = -log10(fdr), color = trend, fill = trend)) +
  ylim(0, max(-log10(results$fdr))+0.5)+
  xlim(-25,25)+
  geom_point(cex = 2.5) + 
  #ggrepel::geom_text_repel(aes(label = tag), size = 2, box.padding = 0.5, max.overlaps = 40,show.legend = FALSE) +
  scale_color_manual(values = alpha(c("Not-sig" = "gray50", "up" = "red", "down" = "blue"), 0.75)) +
  scale_fill_manual(values = c("Not-sig" = "gray50", "up" = "red", "down" = "blue")) +  # 设置填充颜色
  ggrepel::geom_label_repel(aes(label = tag), 
                            size = 2, 
                            box.padding = 0.5, 
                            max.overlaps = 40, 
                            show.legend = FALSE, 
                            color = "white",  # 设置文本颜色为白色
                            label.size = 0,
                            segment.color = "black")+  # 去掉标签边框
  my_theme1 + labs(title = comp)
ggsave(glue(output_dir,"/single_volcano_filtered_padj_mmul.pdf"), vol1, width = 5, height = 5)











## diff sscrofa-----------------------------------------------------------------
file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_sscrofa.csv"
expr_matrix <- read_csv(file_path)
picture_list <- list()


results <- expr_matrix
comp <- "LT_vs_LC_sscrofa"
results <- results %>% as.data.frame() %>% 
  mutate(trend = case_when((PValue < 0.01) & (logFC > 0.5) ~ "up", 
                           (PValue < 0.01) & (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
results <- results[!duplicated(results$gene_name),]
top_diff_genes <- results[order(results$PValue), ][1:10, ] %>% .[.$trend != "Not-sig", "gene_name"]

label_frame <- data_frame(gene_name = c(unique(top_diff_genes)),
                          tag = c(unique(top_diff_genes)))
results <- left_join(results, label_frame, by = "gene_name")
## volcanic
vol1 <- ggplot(results, aes(x = logFC, y = -log10(PValue), color = trend, fill = trend)) +
  ylim(0, max(-log10(results$PValue))+0.5)+
  xlim(-25,25)+
  geom_point(cex = 2.5) + 
  ggrepel::geom_label_repel(aes(label = tag), 
                            size = 2, 
                            box.padding = 0.5, 
                            max.overlaps = 40, 
                            show.legend = FALSE, 
                            color = "white",  # 设置文本颜色为白色
                            label.size = 0,
                            segment.color = "black")+  # 去掉标签边框
  #ggrepel::geom_text_repel(aes(label = tag), size = 2, box.padding = 0.5, max.overlaps = 40,show.legend = FALSE) +
  scale_color_manual(values = alpha(c("Not-sig" = "gray50", "up" = "red", "down" = "blue"), 0.75)) +
  scale_fill_manual(values = c("Not-sig" = "gray50", "up" = "red", "down" = "blue")) +  # 设置填充颜色
  
  my_theme1 + labs(title = comp)


ggsave(glue(output_dir,"/single_volcano_filtered_pval_sscrofa.pdf"), vol1, width = 5, height = 5)


results <- expr_matrix
comp <- "LT_vs_LC_sscrofa"
results <- results %>% as.data.frame() %>% 
  mutate(trend = case_when((fdr < 0.05) & (logFC > 0.5) ~ "up", 
                           (fdr < 0.05) & (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
results <- results[!duplicated(results$gene_name),]
top_diff_genes <- results[order(results$fdr), ][1:10, ] %>% .[.$trend != "Not-sig", "gene_name"]

label_frame <- data_frame(gene_name = c(unique(top_diff_genes)),
                          tag = c(unique(top_diff_genes)))
results <- left_join(results, label_frame, by = "gene_name")
## volcanic
vol1 <- ggplot(results, aes(x = logFC, y = -log10(fdr), color = trend, fill = trend)) +
  ylim(0, max(-log10(results$fdr))+0.5)+
  xlim(-25,25)+
  geom_point(cex = 2.5) + 
  #ggrepel::geom_text_repel(aes(label = tag), size = 2, box.padding = 0.5, max.overlaps = 40,show.legend = FALSE) +
  scale_color_manual(values = alpha(c("Not-sig" = "gray50", "up" = "red", "down" = "blue"), 0.75)) +
  scale_fill_manual(values = c("Not-sig" = "gray50", "up" = "red", "down" = "blue")) +  # 设置填充颜色
  ggrepel::geom_label_repel(aes(label = tag), 
                            size = 2, 
                            box.padding = 0.5, 
                            max.overlaps = 40, 
                            show.legend = FALSE, 
                            color = "white",  # 设置文本颜色为白色
                            label.size = 0,
                            segment.color = "black")+  # 去掉标签边框
  my_theme1 + labs(title = comp)
ggsave(glue(output_dir,"/single_volcano_filtered_padj_sscrofa.pdf"), vol1, width = 5, height = 5)



## diff total-----------------------------------------------------------------

file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_combined.csv"
expr_matrix <- read_csv(file_path)
picture_list <- list()

results <- expr_matrix
names(results) <- c("gene_name","logFC","logCPM","PValue","fdr")
comp <- "LT_vs_LC_total"
results <- results %>% as.data.frame() %>% 
  mutate(trend = case_when((PValue < 0.01) & (logFC > 0.5) ~ "up", 
                           (PValue < 0.01) & (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
results <- results[!duplicated(results$gene_name),]
top_diff_genes <- results[order(results$PValue), ][1:10, ] %>% .[.$trend != "Not-sig", "gene_name"]

label_frame <- data_frame(gene_name = c(unique(top_diff_genes)),
                          tag = c(unique(top_diff_genes)))
results <- left_join(results, label_frame, by = "gene_name")
## volcanic
vol1 <- ggplot(results, aes(x = logFC, y = -log10(PValue), color = trend, fill = trend)) +
  ylim(0, max(-log10(results$PValue))+0.5)+
  xlim(-25,25)+
  geom_point(cex = 2.5) + 
  ggrepel::geom_label_repel(aes(label = tag), 
                            size = 2, 
                            box.padding = 0.5, 
                            max.overlaps = 40, 
                            show.legend = FALSE, 
                            color = "white",  # 设置文本颜色为白色
                            label.size = 0,
                            segment.color = "black")+  # 去掉标签边框
  #ggrepel::geom_text_repel(aes(label = tag), size = 2, box.padding = 0.5, max.overlaps = 40,show.legend = FALSE) +
  scale_color_manual(values = alpha(c("Not-sig" = "gray50", "up" = "red", "down" = "blue"), 0.75)) +
  scale_fill_manual(values = c("Not-sig" = "gray50", "up" = "red", "down" = "blue")) +  # 设置填充颜色
  
  my_theme1 + labs(title = comp)


ggsave(glue(output_dir,"/single_volcano_filtered_pval_total.pdf"), vol1, width = 5, height = 5)


results <- expr_matrix
names(results) <- c("gene_name","logFC","logCPM","PValue","fdr")
comp <- "LT_vs_LC_total"
results <- results %>% as.data.frame() %>% 
  mutate(trend = case_when((fdr < 0.05) & (logFC > 0.5) ~ "up", 
                           (fdr < 0.05) & (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
results <- results[!duplicated(results$gene_name),]
top_diff_genes <- results[order(results$fdr), ][1:10, ] %>% .[.$trend != "Not-sig", "gene_name"]

label_frame <- data_frame(gene_name = c(unique(top_diff_genes)),
                          tag = c(unique(top_diff_genes)))
results <- left_join(results, label_frame, by = "gene_name")
## volcanic
vol1 <- ggplot(results, aes(x = logFC, y = -log10(fdr), color = trend, fill = trend)) +
  ylim(0, max(-log10(results$fdr))+0.5)+
  xlim(-25,25)+
  geom_point(cex = 2.5) + 
  #ggrepel::geom_text_repel(aes(label = tag), size = 2, box.padding = 0.5, max.overlaps = 40,show.legend = FALSE) +
  scale_color_manual(values = alpha(c("Not-sig" = "gray50", "up" = "red", "down" = "blue"), 0.75)) +
  scale_fill_manual(values = c("Not-sig" = "gray50", "up" = "red", "down" = "blue")) +  # 设置填充颜色
  ggrepel::geom_label_repel(aes(label = tag), 
                            size = 2, 
                            box.padding = 0.5, 
                            max.overlaps = 40, 
                            show.legend = FALSE, 
                            color = "white",  # 设置文本颜色为白色
                            label.size = 0,
                            segment.color = "black")+  # 去掉标签边框
  my_theme1 + labs(title = comp)
ggsave(glue(output_dir,"/single_volcano_filtered_padj_total.pdf"), vol1, width = 5, height = 5)
