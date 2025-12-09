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




##MMUL-----------------------------------------------------------------------------
count_total <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm_mmul.txt")

# 去掉有0的行
count_total <- count_total %>% filter(across(everything(), ~ . != 0))
# 去掉全是 0 的行
#count_total <- count_total[rowSums(count_total != 0) > 0, ]

normalized_counts <- log2(count_total+1)
normalized_counts <- as.data.frame(normalized_counts)

#foldchange
file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_mmul.csv"
fc_matrix <- read_csv(file_path)
genes <- fc_matrix$gene_name
normalized_counts <- normalized_counts[rownames(normalized_counts) %in% genes,]
fc_matrix <- fc_matrix[fc_matrix$gene_name %in% rownames(normalized_counts),]
fc_matrix <- fc_matrix[order(fc_matrix$logFC, decreasing = TRUE), ]
normalized_counts <- normalized_counts[fc_matrix$gene_name,]
fc <- fc_matrix$logFC
names(fc) <- rownames(fc_matrix)
# 定义fold change的颜色映射
fc_col_fun <- colorRamp2(c(min(fc), 0, max(fc)), c("blue", "white", "red"))

# 创建行注释对象（在右侧显示）
ha <- rowAnnotation(
  FoldChange = fc,  # 这里传入fc向量，anno_simple会自动进行颜色映射
  col = list(FoldChange = fc_col_fun),
  annotation_name_side = "top",  # 注释名称放在顶部
  annotation_legend_param = list(title = "log2FC")
)
count_qual <- c(normalized_counts$LC,normalized_counts$LT)
quintiles <- quantile(count_qual, probs = c(0, 0.2, 0.5, 0.8, 1))
print(quintiles)

#heatmap
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/total_genes_mmul_heatmap_v2_1.pdf", width = 3, height = 5)
Heatmap(normalized_counts,
        name = "Expression_mmul",
        #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
        col = colorRamp2(c(0, 0.15, 1), c("#3B035F", "#07989E", "#F8D33E")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        column_title = NULL,
        right_annotation = ha,  # 将行注释放在右侧
        column_split = 9)
dev.off()

## all zero-------------------------------------------------------------------
count_total <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm_mmul.txt")

# 去掉有0的行
#count_total <- count_total %>% filter(across(everything(), ~ . != 0))
# 去掉全是 0 的行
count_total <- count_total[rowSums(count_total != 0) > 0, ]

#normalized_counts <- count_total * 1e6
normalized_counts <- log2(count_total+1)
normalized_counts <- as.data.frame(normalized_counts)

#foldchange
file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_mmul.csv"
fc_matrix <- read_csv(file_path)
genes <- fc_matrix$gene_name
normalized_counts <- normalized_counts[rownames(normalized_counts) %in% genes,]
fc_matrix <- fc_matrix[fc_matrix$gene_name %in% rownames(normalized_counts),]
fc_matrix <- fc_matrix[order(fc_matrix$logFC, decreasing = TRUE), ]
normalized_counts <- normalized_counts[fc_matrix$gene_name,]
fc <- fc_matrix$logFC
names(fc) <- rownames(fc_matrix)
# 定义fold change的颜色映射
fc_col_fun <- colorRamp2(c(min(fc), 0, max(fc)), c("blue", "white", "red"))

# 创建行注释对象（在右侧显示）
ha <- rowAnnotation(
  FoldChange = fc,  # 这里传入fc向量，anno_simple会自动进行颜色映射
  col = list(FoldChange = fc_col_fun),
  annotation_name_side = "top",  # 注释名称放在顶部
  annotation_legend_param = list(title = "log2FC")
)


#total
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/total_genes_mmul_heatmap_v2_1_all_zero_exclude.pdf", width = 3, height = 5)
Heatmap(normalized_counts,
        name = "Expression_mmul",
        #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
        col = colorRamp2(c(0, 0.15, 1), c("#3B035F", "#07989E", "#F8D33E")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        column_title = NULL,
        right_annotation = ha,  # 将行注释放在右侧
        column_split = 9)
dev.off()













##sscrofa-----------------------------------------------------------------------------
count_total <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm_sscrofa.txt")
# 去掉有0的行
count_total <- count_total %>% filter(across(everything(), ~ . != 0))
# 去掉全是 0 的行
#count_total <- count_total[rowSums(count_total != 0) > 0, ]

normalized_counts <- log2(normalized_counts+1)
normalized_counts <- as.data.frame(normalized_counts)

#foldchange
file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_sscrofa.csv"
fc_matrix <- read_csv(file_path)

genes <- fc_matrix$gene_name
normalized_counts <- normalized_counts[rownames(normalized_counts) %in% genes,]
fc_matrix <- fc_matrix[fc_matrix$gene_name %in% rownames(normalized_counts),]
fc_matrix <- fc_matrix[order(fc_matrix$logFC, decreasing = TRUE), ]
normalized_counts <- normalized_counts[fc_matrix$gene_name,]
fc <- fc_matrix$logFC
names(fc) <- fc_matrix$gene_name
# 定义fold change的颜色映射
fc_col_fun <- colorRamp2(c(min(fc), 0, max(fc)), c("blue", "white", "red"))

# 创建行注释对象（在右侧显示）
ha <- rowAnnotation(
  FoldChange = fc,  # 这里传入fc向量，anno_simple会自动进行颜色映射
  col = list(FoldChange = fc_col_fun),
  annotation_name_side = "top",  # 注释名称放在顶部
  annotation_legend_param = list(title = "log2FC")
)

count_qual <- c(normalized_counts$LC,normalized_counts$LT)
quintiles <- quantile(count_qual, probs = c(0, 0.2, 0.5, 0.8, 1))
print(quintiles)



pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/total_genes_sscrofa_heatmap_v2_1.pdf", width = 3, height = 5)
Heatmap(normalized_counts,
        name = "Expression_sscrofa",
        #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
        col = colorRamp2(c(0, 0.15, 1), c("#3B035F", "#07989E", "#F8D33E")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        column_title = NULL,
        right_annotation = ha,  # 将行注释放在右侧
        column_split = 9)
dev.off()


##all_zero_exclude-------------------------------------------------------------
count_total <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm_sscrofa.txt")
# 去掉有0的行
#count_total <- count_total %>% filter(across(everything(), ~ . != 0))
# 去掉全是 0 的行
count_total <- count_total[rowSums(count_total != 0) > 0, ]

normalized_counts <- count_total * 1e6
normalized_counts <- log10(normalized_counts+1)
normalized_counts <- as.data.frame(normalized_counts)

#foldchange
file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_sscrofa.csv"
fc_matrix <- read_csv(file_path)

genes <- fc_matrix$gene_name
normalized_counts <- normalized_counts[rownames(normalized_counts) %in% genes,]
fc_matrix <- fc_matrix[fc_matrix$gene_name %in% rownames(normalized_counts),]
fc_matrix <- fc_matrix[order(fc_matrix$logFC, decreasing = TRUE), ]
normalized_counts <- normalized_counts[fc_matrix$gene_name,]
fc <- fc_matrix$logFC
names(fc) <- fc_matrix$gene_name
# 定义fold change的颜色映射
fc_col_fun <- colorRamp2(c(min(fc), 0, max(fc)), c("blue", "white", "red"))

# 创建行注释对象（在右侧显示）
ha <- rowAnnotation(
  FoldChange = fc,  # 这里传入fc向量，anno_simple会自动进行颜色映射
  col = list(FoldChange = fc_col_fun),
  annotation_name_side = "top",  # 注释名称放在顶部
  annotation_legend_param = list(title = "log2FC")
)




pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/total_genes_sscrofa_heatmap_v2_1_all_zero_exclude.pdf", width = 3, height = 5)
Heatmap(normalized_counts,
        name = "Expression_sscrofa",
        #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
        col = colorRamp2(c(0, 4, 8), c("#3B035F", "#07989E", "#F8D33E")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        column_title = NULL,
        right_annotation = ha,  # 将行注释放在右侧
        column_split = 9)
dev.off()
