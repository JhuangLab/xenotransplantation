.libPaths(new = c(.libPaths(), "~/sbin/R/R-4.3.0", "~/sbin/R/R-4.2.1"))
# tryCatch({library(clusterProfiler)}, error = \(e) library(clusterProfiler))
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "QFeatures", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "MSnbase", "msqrob2", 
          "enrichplot", "org.Hs.eg.db", "parallel", "jhuanglabRNAseq", "jhuanglabGO", "limma", 
          "readxl", "writexl", "ComplexHeatmap", "circlize", "clusterProfiler", "MsCoreUtils", 
          "missForest", "paletteer", "factoextra", "FactoMineR","VennDiagram")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "liver"
dataset <- "wangrongrong"
species <- "Xeno_pig_monkey"
workdir <- glue::glue("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn")
workdir %>% fs::dir_create() %>% setwd()


up_regulate_list <- list()
down_regulate_list <- list()

##mmul regulate genes
file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_mmul.csv"
expr_matrix_mmul <- read_csv(file_path)
expr_matrix_mmul <- expr_matrix_mmul %>% as.data.frame() %>% 
  mutate(trend = case_when((fdr < 0.05) & (logFC > 0.5) ~ "up", 
                           (fdr < 0.05) & (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
up_regulate_list[["up_regulate_mmul"]] <- expr_matrix_mmul[expr_matrix_mmul$trend == "up", c("gene_name")] %>% na.omit()
down_regulate_list[["down_regulate_mmul"]] <- expr_matrix_mmul[expr_matrix_mmul$trend == "down", c("gene_name")] %>% na.omit()

##sscrofa regulate genes
file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_sscrofa.csv"
expr_matrix_sscrofa <- read_csv(file_path)
expr_matrix_sscrofa <- expr_matrix_sscrofa %>% as.data.frame() %>% 
  mutate(trend = case_when((fdr < 0.05) & (logFC > 0.5) ~ "up", 
                           (fdr < 0.05) & (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
up_regulate_list[["up_regulate_sscrofa"]] <- expr_matrix_sscrofa[expr_matrix_sscrofa$trend == "up", c("gene_name")] %>% na.omit()
down_regulate_list[["down_regulate_sscrofa"]] <- expr_matrix_sscrofa[expr_matrix_sscrofa$trend == "down", c("gene_name")] %>% na.omit()

fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
expr_matrix_protein <- fc_table %>% as.data.frame() %>% 
  mutate(trend = case_when((FDR < 0.05) & (logFC > 0.5) ~ "up", 
                           (FDR < 0.05) & (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
up_regulate_list[["up_regulate_protein"]] <- expr_matrix_protein[expr_matrix_protein$trend == "up", ] %>% rownames() %>% na.omit()
down_regulate_list[["down_regulate_protein"]]  <- expr_matrix_protein[expr_matrix_protein$trend == "down", ] %>% rownames() %>% na.omit()

intersection_ABC_up <- Reduce(intersect, up_regulate_list)  # A、B 和 C 的交集
write.csv(intersection_ABC_up, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/up_all.csv")

intersect1 <- intersect(up_regulate_list[["up_regulate_mmul"]], up_regulate_list[["up_regulate_protein"]])
write.csv(intersect1, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/up_mmul_protein.csv")
intersect2 <- intersect(up_regulate_list[["up_regulate_sscrofa"]], up_regulate_list[["up_regulate_protein"]])
write.csv(intersect2, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/up_sscrofa_protein.csv")
intersect3 <- intersect(up_regulate_list[["up_regulate_mmul"]], up_regulate_list[["up_regulate_sscrofa"]])
write.csv(intersect3, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/up_mmul_sscrofa.csv")

intersection_ABC_down <- Reduce(intersect, down_regulate_list)  # A、B 和 C 的交集
write.csv(intersection_ABC_down, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/down_all.csv")

intersect1 <- intersect(down_regulate_list[["down_regulate_mmul"]], down_regulate_list[["down_regulate_protein"]])
write.csv(intersect1, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/down_mmul_protein.csv")
intersect2 <- intersect(down_regulate_list[["down_regulate_sscrofa"]], down_regulate_list[["down_regulate_protein"]])
write.csv(intersect2, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/down_sscrofa_protein.csv")
intersect3 <- intersect(down_regulate_list[["down_regulate_mmul"]], down_regulate_list[["down_regulate_sscrofa"]])
write.csv(intersect3, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/down_mmul_sscrofa.csv")


venn.plot <- venn.diagram(
  x = up_regulate_list,
  category.names = names(up_regulate_list),
  filename = NULL,  # 设置为NULL以在R图形设备中绘制
  output = TRUE,
  alpha = 0.5,  # 设置透明度
  cat.cex = 1,  # 设置类别标签的字体大小
  fill = c("#D72528", "#2AA12B", "#FD7E0D"),  # 设置圆圈颜色
  cex = 1,  # 设置数字的字体大小
  cat.pos = c(335, 25, 180),  # 设置类别标签的位置，单位为度
  cat.dist = c(0.06, 0.06, 0.05)# 设置类别标签与圆圈的距离
)
ggsave(venn.plot,file="/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/venn_up_regulate.pdf", width = 6 ,height = 6)

venn.plot <- venn.diagram(
  x = down_regulate_list,
  category.names = names(down_regulate_list),
  filename = NULL,  # 设置为NULL以在R图形设备中绘制
  output = TRUE,
  alpha = 0.5,  # 设置透明度
  cat.cex = 1,  # 设置类别标签的字体大小
  fill = c("#D72528", "#2AA12B", "#FD7E0D"),  # 设置圆圈颜色
  cex = 1,  # 设置数字的字体大小
  cat.pos = c(335, 25, 180),  # 设置类别标签的位置，单位为度
  cat.dist = c(0.06, 0.06, 0.05)# 设置类别标签与圆圈的距离
)
ggsave(venn.plot,file="/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/venn_down_regulate.pdf", width = 6 ,height = 6)








# only FC
up_regulate_list <- list()
down_regulate_list <- list()

##mmul regulate genes
file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_mmul.csv"
expr_matrix_mmul <- read_csv(file_path)
expr_matrix_mmul <- expr_matrix_mmul %>% as.data.frame() %>% 
  mutate(trend = case_when((logFC > 0.5) ~ "up", 
                           (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
up_regulate_list[["up_regulate_mmul"]] <- expr_matrix_mmul[expr_matrix_mmul$trend == "up", c("gene_name")] %>% na.omit()
down_regulate_list[["down_regulate_mmul"]] <- expr_matrix_mmul[expr_matrix_mmul$trend == "down", c("gene_name")] %>% na.omit()

##sscrofa regulate genes
file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_sscrofa.csv"
expr_matrix_sscrofa <- read_csv(file_path)
expr_matrix_sscrofa <- expr_matrix_sscrofa %>% as.data.frame() %>% 
  mutate(trend = case_when((logFC > 0.5) ~ "up", 
                           (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
up_regulate_list[["up_regulate_sscrofa"]] <- expr_matrix_sscrofa[expr_matrix_sscrofa$trend == "up", c("gene_name")] %>% na.omit()
down_regulate_list[["down_regulate_sscrofa"]] <- expr_matrix_sscrofa[expr_matrix_sscrofa$trend == "down", c("gene_name")] %>% na.omit()

fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
expr_matrix_protein <- fc_table %>% as.data.frame() %>% 
  mutate(trend = case_when((logFC > 0.5) ~ "up", 
                           (logFC < -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
up_regulate_list[["up_regulate_protein"]] <- expr_matrix_protein[expr_matrix_protein$trend == "up", ] %>% rownames() %>% na.omit()
down_regulate_list[["down_regulate_protein"]]  <- expr_matrix_protein[expr_matrix_protein$trend == "down", ] %>% rownames() %>% na.omit()

intersection_ABC_up <- Reduce(intersect, up_regulate_list)  # A、B 和 C 的交集
intersection_ABC_down <- Reduce(intersect, down_regulate_list)  # A、B 和 C 的交集
# 计算交集
intersection_AB <- intersect(up_regulate_list[["up_regulate_mmul"]], up_regulate_list[["up_regulate_sscrofa"]])  # A 和 B 的交集
intersection_AC <- intersect(up_regulate_list[["up_regulate_mmul"]], up_regulate_list[["up_regulate_protein"]])  # A 和 C 的交集
intersection_BC <- intersect(up_regulate_list[["up_regulate_sscrofa"]], up_regulate_list[["up_regulate_protein"]])  # B 和 C 的交集
intersection_ABC <- Reduce(intersect, list(setA, setB, setC))  # A、B 和 C 的交集

# 创建一个数据框以存储交集结果
gene_intersections <- data.frame(
  Gene = c(intersection_AB, intersection_AC, intersection_BC, intersection_ABC),
  Source = c(rep("A ∩ B", length(intersection_AB)),
             rep("A ∩ C", length(intersection_AC)),
             rep("B ∩ C", length(intersection_BC)),
             rep("A ∩ B ∩ C", length(intersection_ABC)))
)

# 输出交集的基因列表
print(gene_intersections)



venn.plot <- venn.diagram(
  x = up_regulate_list,
  category.names = names(up_regulate_list),
  filename = NULL,  # 设置为NULL以在R图形设备中绘制
  output = TRUE,
  alpha = 0.5,  # 设置透明度
  cat.cex = 1,  # 设置类别标签的字体大小
  fill = c("#D72528", "#2AA12B", "#FD7E0D"),  # 设置圆圈颜色
  cex = 1,
  height = 5,  # 设置高度
  width = 5,   # 设置宽度
  # 设置数字的字体大小
  cat.pos = c(335, 25, 180),  # 设置类别标签的位置，单位为度
  cat.dist = c(0.06, 0.06, 0.05)# 设置类别标签与圆圈的距离
)
ggsave(venn.plot,file="/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/venn_up_regulate_v2_test.pdf", width = 6 ,height = 6)

venn.plot <- venn.diagram(
  x = down_regulate_list,
  category.names = names(down_regulate_list),
  filename = NULL,  # 设置为NULL以在R图形设备中绘制
  output = TRUE,
  alpha = 0.5,  # 设置透明度
  cat.cex = 1,  # 设置类别标签的字体大小
  fill = c("#D72528", "#2AA12B", "#FD7E0D"),  # 设置圆圈颜色
  cex = 1,  # 设置数字的字体大小
  cat.pos = c(335, 25, 180),  # 设置类别标签的位置，单位为度
  cat.dist = c(0.06, 0.06, 0.05)# 设置类别标签与圆圈的距离
)
ggsave(venn.plot,file="/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/venn_down_regulate_v2.pdf", width = 6 ,height = 6)
