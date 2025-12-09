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

all_exp <- read.delim("/cluster/home/jhuang/projects/liver/data/zhangwei/pig/protein/report/report-liver/data/3.comparison/PIT_CON-IT/PIT_CON-IT.all.xls")
all_exp <- all_exp[,c("gene_name","PIT","CON.IT")]
all_exp <- all_exp[!grepl("-", all_exp$gene_name), ] %>% na.omit()
all_exp <- all_exp[!duplicated(all_exp$gene_name),]
rownames(all_exp) <- NULL
all_exp <- all_exp %>% column_to_rownames(.,"gene_name")

column_sums <- colSums(all_exp)
normalized_counts <- t(t(all_exp) / column_sums) * 1e6  # 将标准化后的值乘以 1,000,000
normalized_counts <- log10(normalized_counts+1)
normalized_counts <- as.data.frame(normalized_counts)
fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)

# 获取higher in treatment group

sorted_indices <- order(fc_table[, "logFC"], decreasing = TRUE)
top_indices <- sorted_indices[1:25]
# 保留这些索引所在的行
treat_rows <- fc_table[top_indices, ] %>% rownames()

##获取 lower in ctrl group
sorted_indices <- order(fc_table[, "logFC"])
top_indices <- sorted_indices[1:25]
# 保留这些索引所在的行
ctrl_rows <- fc_table[top_indices, ] %>% rownames()

normalized_counts <- normalized_counts[c(treat_rows,ctrl_rows),]
# 去掉特定的行
filtered_counts <- normalized_counts[!rownames(normalized_counts) %in% c("FGA", "FGB", "FGG", "HBB"), ]
# 标准化数据
#scaled_data <- t(apply(pr_dat, 1, scale)) %>% as.data.frame()  # 按行标准化
#names(scaled_data) <- names(pr_dat)


pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/all_gene_heatmap_25_2.pdf", width = 4, height = 5)
Heatmap(filtered_counts,
        name = "Expression_all",
        col = colorRampPalette(c("#3B035F", "#07989E", "#F8D33E"))(50),
        #col = colorRamp2(c(0, 1, 2), c("#3B035F", "#07989E", "#F8D33E")),
        show_row_names = TRUE,
        show_column_names = TRUE,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        column_title = NULL,
        row_names_gp = gpar(fontsize = 7), # 设置行名字体大小为 10
        column_split = 9)
dev.off()


##v2_2
all_exp <- read.delim("/cluster/home/jhuang/projects/liver/data/zhangwei/pig/protein/report/report-liver/data/3.comparison/PIT_CON-IT/PIT_CON-IT.all.xls")
gene_names <- c("ITGB1", "ITGA5", "ITGB3", "FN1", "PECAM1","CALR","MAP2K1",
                "IDH1","PSAT1","MIF")
all_exp <- all_exp[,c("gene_name","PIT","CON.IT")]
all_exp <- all_exp[!grepl("-", all_exp$gene_name), ] %>% na.omit()
all_exp <- all_exp[!duplicated(all_exp$gene_name),]
rownames(all_exp) <- NULL
all_exp <- all_exp %>% column_to_rownames(.,"gene_name")

column_sums <- colSums(all_exp)
normalized_counts <- t(t(all_exp) / column_sums) * 1e6  # 将标准化后的值乘以 1,000,000
normalized_counts <- log10(normalized_counts+1)
normalized_counts <- as.data.frame(normalized_counts)
normalized_counts <- normalized_counts %>% filter(across(everything(), ~ . != 0))

fc_matrix <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
fc_matrix$gene_name <- rownames(fc_matrix)
# 获取higher in treatment group
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
count_qual <- c(normalized_counts$PIT,normalized_counts$CON.IT)
quintiles <- quantile(count_qual, probs = c(0, 0.2, 0.5, 0.8, 1))
print(quintiles)

ht <- Heatmap(normalized_counts,
              name = "Expression_protein",
              #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
              col = colorRamp2(c(0, 1.25, 2.5), c("#3B035F", "#07989E", "#F8D33E")),
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
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/all_gene_heatmap_v2_2_TEST.pdf", width = 3.5, height = 5)
draw(ht)
dev.off()



##v2 all gene 0 exclude
all_exp <- read.delim("/cluster/home/jhuang/projects/liver/data/zhangwei/pig/protein/report/report-liver/data/3.comparison/PIT_CON-IT/PIT_CON-IT.all.xls")
all_exp <- all_exp[,c("gene_name","PIT","CON.IT")]
all_exp <- all_exp[!grepl("-", all_exp$gene_name), ] %>% na.omit()
all_exp <- all_exp[!duplicated(all_exp$gene_name),]
rownames(all_exp) <- NULL
all_exp <- all_exp %>% column_to_rownames(.,"gene_name")
fc_matrix$gene_name <- rownames(fc_matrix)

column_sums <- colSums(all_exp)
normalized_counts <- t(t(all_exp) / column_sums) * 1e6  # 将标准化后的值乘以 1,000,000
normalized_counts <- log10(normalized_counts+1)
normalized_counts <- as.data.frame(normalized_counts)

count_total <- count_total[rowSums(count_total != 0) > 0, ]

fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)

# 获取higher in treatment group
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


pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/all_gene_heatmap_v2_1_all_zero_exclude.pdf", width = 3, height = 5)
Heatmap(normalized_counts,
        name = "Expression_protein",
        #col = colorRampPalette(c("#003366", "#FFFFFF","#990000"))(200),
        col = colorRamp2(c(0, 1.25, 2.5), c("#3B035F", "#07989E", "#F8D33E")),
        show_row_names = FALSE,
        show_column_names = TRUE,
        cluster_rows = FALSE,
        cluster_columns = TRUE,
        column_title = NULL,
        right_annotation = ha,  # 将行注释放在右侧
        column_split = 9)
dev.off()





