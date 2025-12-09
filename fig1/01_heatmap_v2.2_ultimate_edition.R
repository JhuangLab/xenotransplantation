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


##total-------------------------------------------------------------------------
count_total <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm.txt")
# 去掉全是 0 的行
count_total <- count_total[rowSums(count_total != 0) > 0, ]

#high variance gene 2000
var_phos <- apply(count_total, 1, var)
phosps <- sort(-var_phos)[1:2000] %>% names() %>% unique()
fpkm_fn_heatmap_2000 <- count_total[phosps, ]
log_frame <- data.frame(
  LC = log10(fpkm_fn_heatmap_2000$LC + 1),
  LT = log10(fpkm_fn_heatmap_2000$LT + 1)
)
rownames(log_frame) <- rownames(fpkm_fn_heatmap_2000)

pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/2000_genes_all.pdf", width = 3, height = 5)
Heatmap(log_frame,
        name = "Expression_all",
        #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
        col = colorRamp2(c(0, 1, 2), c("#3B035F", "#07989E", "#F8D33E")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = NULL,
        column_split = 9)
dev.off()



##MMUL-----------------------------------------------------------------------------
count_total <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm_mmul.txt")
gene_names <- c("VCAN", "FN1", "ARF1", "CTSS", "CD74","TGFBI",
                "SLIT2","TEAD1","MAST4")# Replace with your actual gene names
# 去掉有0的行
count_total <- count_total %>% filter(across(everything(), ~ . != 0))
# 去掉全是 0 的行
#count_total <- count_total[rowSums(count_total != 0) > 0, ]

normalized_counts <- count_total * 1e6
normalized_counts <- log10(normalized_counts+1)
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

ht <- Heatmap(normalized_counts,
              name = "Expression_mmul",
              #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
              col = colorRamp2(c(0, 3, 5), c("#3B035F", "#07989E", "#F8D33E")),
              show_row_names = FALSE,
              show_column_names = TRUE,
              cluster_rows = FALSE,
              cluster_columns = TRUE,
              column_title = NULL,
              right_annotation = ha,  # 将行注释放在右侧
              column_split = 9)



idx_gn <- which(rownames(normalized_counts) %in% gene_names)
gn_lab <- rownames(normalized_counts)[idx_gn]
row_annotation <- rowAnnotation(link = anno_mark(at = idx_gn, 
                                                 labels = gn_lab, labels_gp = gpar(fontsize = 10),
                                                 side = "left"))
ht <- row_annotation + ht

#heatmap
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/total_genes_mmul_heatmap_v2_2.pdf", width = 3.5, height = 5)
draw(ht)
dev.off()

## all zero-------------------------------------------------------------------
count_total <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm_mmul.txt")

# 去掉有0的行
#count_total <- count_total %>% filter(across(everything(), ~ . != 0))
# 去掉全是 0 的行
count_total <- count_total[rowSums(count_total != 0) > 0, ]

normalized_counts <- count_total * 1e6
normalized_counts <- log10(normalized_counts+1)
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
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/total_genes_mmul_heatmap_v2_all_zero_exclude.pdf", width = 3, height = 5)
Heatmap(normalized_counts,
        name = "Expression_mmul",
        #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
        col = colorRamp2(c(0, 3, 5), c("#3B035F", "#07989E", "#F8D33E")),
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
gene_names <- c("CXCL14", "CCL5", "TNFSF13B", "XCR1","ITGA4", "CD4", "CD8A","CD74","MCM3","MCM4","SFN","SGO1","COQ7")

# 去掉有0的行
count_total <- count_total %>% filter(across(everything(), ~ . != 0))
# 去掉全是 0 的行
#count_total <- count_total[rowSums(count_total != 0) > 0, ]

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

ht <- Heatmap(normalized_counts,
        name = "Expression_sscrofa",
        #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
        col = colorRamp2(c(0, 4, 8), c("#3B035F", "#07989E", "#F8D33E")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = F,
        cluster_columns = TRUE,
        column_title = NULL,
        right_annotation = ha,
        #left_annotation = gene_annotation,
        # 将行注释放在右侧
        column_split = 9)
idx_gn <- which(rownames(normalized_counts) %in% gene_names)
gn_lab <- rownames(normalized_counts)[idx_gn]
row_annotation <- rowAnnotation(link = anno_mark(at = idx_gn, 
                                   labels = gn_lab, labels_gp = gpar(fontsize = 10),
                                   side = "left"))
ht <- row_annotation + ht

pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/total_genes_sscrofa_heatmap_v2_2.pdf", width = 4, height = 5)
draw(ht)
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




pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/total_genes_sscrofa_heatmap_v2_all_zero_exclude.pdf", width = 3, height = 5)
Heatmap(normalized_counts,
        name = "Expression_sscrofa",
        #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
        col = colorRamp2(c(0, 4, 8), c("#3B035F", "#07989E", "#F8D33E")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        column_title = NULL,
        right_annotation = ha,
        left_annotation = mark_annotation,# 将行注释放在右侧
        column_split = 9)
dev.off()
