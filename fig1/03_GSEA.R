.libPaths(new = c(.libPaths(), "/cluster/home/danyang_jh/sbin/R/R-4.3.0", "/cluster/home/danyang_jh/sbin/R/R-4.2.1","/cluster/home/xyzhang_jh/sbin/R"))
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "QFeatures", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "MSnbase", "msqrob2", 
          "enrichplot", "org.Hs.eg.db", "parallel", "jhuanglabRNAseq", "jhuanglabGO", "limma", 
          "readxl", "writexl", "ComplexHeatmap", "circlize", "clusterProfiler", "MsCoreUtils", 
          "missForest", "paletteer", "factoextra", "FactoMineR", "openxlsx","gdata","ggrepel","clusterProfiler","org.Ss.eg.db","DOSE","org.Mmu.eg.db")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "liver"
dataset <- "wangrongrong"
species <- "Xeno_pig_monkey"
workdir <- glue::glue("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA")
workdir %>% fs::dir_create() %>% setwd()



## mmul -------------------------------------------------------------------------
df <- read.csv("/cluster/home/xyzhang_jh/projects/Xenograft/monkey/diff_genes_mmul.csv")
df <- df[,c("gene_name", "logFC", "PValue")]



# 将基因名转换为 Entrez ID
gene_entrez <- bitr(df$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mmu.eg.db")
df_mmul_sig <- merge(df, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL")

# 过滤掉没有对应 Entrez ID 的基因
df_mmul_sig <- na.omit(df_mmul_sig)

# 创建一个基因排名向量
gene_list <- df_mmul_sig$logFC
names(gene_list) <- df_mmul_sig$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
saveRDS(gene_list, "GSEA_mmu_gene_list.rds")

#fat1
gene_list <- read_rds(file = "/cluster/home/xyzhang_jh/projects/Xenograft/monkey/GSEA_mmu_gene_list.rds")
gene_list <- read_rds(file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/monkey/GSEA_mmu_gene_list.rds")

# 运行 GSEA
gsea_results <- gseKEGG(geneList = gene_list,
                        organism = "mcc",  # 人类
                        nPerm = 1000,      # 随机化次数
                        minGSSize = 10,    # 基因集最小大小
                        maxGSSize = 500,   # 基因集最大大小
                        pvalueCutoff = 1,
                        verbose = TRUE)

saveRDS(gsea_results,"/cluster/home/xyzhang_jh/projects/Xenograft/monkey/GSEA_monkey_kegg.rds")


#node70
gsea_results_kegg <- readRDS(file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/monkey/GSEA_monkey_kegg.rds")
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/monkey/facet_by_sign_kegg.pdf",width = 8,height = 10)
dotplot(gsea_results_kegg,split=".sign")+facet_grid(~.sign)+
  theme(
    axis.text = element_text(size = 8),  # 坐标轴文本大小
    axis.title = element_text(size = 12), # 坐标轴标题大小
    strip.text = element_text(size = 10), # 分面标签大小
    plot.title = element_text(size = 14)  # 图表标题大小
  ) +
  ggtitle("GSEA KEGG Analysis LT_vs_LC")  # 添加标题
dev.off()

gsea_results_kegg2 <- as.data.frame(gsea_results_kegg)
convert_genes <- function(gene_ids) {
  ids <- unlist(strsplit(gene_ids, "/"))
  gene_info <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Mmu.eg.db")
}
# 应用转换函数到 core_enrichment 列
a <- sapply(gsea_results_kegg2$core_enrichment, convert_genes)

gsea_results_kegg2$core_enrichment <- sapply(a[2,], function(x) {
  paste(unlist(x), collapse = "/")
})
write.csv(gsea_results_kegg2, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/monkey/gsea_kegg_res.csv")

# 运行 GSEA
gsea_results_go <- gseGO(geneList = gene_list,
                      OrgDb = "org.Mmu.eg.db",  # 人类
                      nPerm = 1000,      # 随机化次数
                      minGSSize = 10,    # 基因集最小大小
                      maxGSSize = 500,   # 基因集最大大小
                      pvalueCutoff = 1,
                      verbose = TRUE)
saveRDS(gsea_results_go,"/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/monkey/GSEA_monkey_go.rds")
gsea_results_go <- readRDS("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/monkey/GSEA_monkey_go.rds")
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/monkey/facet_by_sign_go.pdf",width = 8,height = 10)
dotplot(gsea_results_go,split=".sign")+facet_grid(~.sign)+
  theme(
    axis.text = element_text(size = 8),  # 坐标轴文本大小
    axis.title = element_text(size = 12), # 坐标轴标题大小
    strip.text = element_text(size = 10), # 分面标签大小
    plot.title = element_text(size = 14)  # 图表标题大小
  ) +
  ggtitle("GSEA GO Analysis LT_vs_LC")  # 添加标题
dev.off()

gsea_results_go2 <- as.data.frame(gsea_results_go)
convert_genes <- function(gene_ids) {
  ids <- unlist(strsplit(gene_ids, "/"))
  gene_info <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Mmu.eg.db")
}
# 应用转换函数到 core_enrichment 列
a <- sapply(gsea_results_go2$core_enrichment, convert_genes)

gsea_results_go2$core_enrichment <- sapply(a[2,], function(x) {
  paste(unlist(x), collapse = "/")
})
write.csv(gsea_results_go2, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/monkey/gsea_go_res.csv")




#sscrofa------------------------------------------------------------------------
df <- read.csv("/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/diff/diff_genes_sscrofa.csv")
df <- df[,c("gene_name", "logFC", "PValue")]
# 示例数据框，包含基因名和对应的统计量
# 假设 df 是您的差异表达结果数据框


# 将基因名转换为 Entrez ID
gene_entrez <- bitr(df$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Ss.eg.db")
df_mmul_sig <- merge(df, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL")

# 过滤掉没有对应 Entrez ID 的基因
df_mmul_sig <- na.omit(df_mmul_sig)

# 创建一个基因排名向量
gene_list <- df_mmul_sig$logFC
names(gene_list) <- df_mmul_sig$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)
saveRDS(gene_list, "GSEA_sscrofa_gene_list.rds")

#fat1
gene_list <- read_rds(file = "/cluster/home/xyzhang_jh/projects/Xenograft/pig/GSEA_sscrofa_gene_list.rds")
# 运行 GSEA
gsea_results <- gseKEGG(geneList = gene_list,
                        organism = "ssc",  # 猴
                        nPerm = 1000,      # 随机化次数
                        minGSSize = 10,    # 基因集最小大小
                        maxGSSize = 500,   # 基因集最大大小
                        pvalueCutoff = 1,
                        verbose = TRUE)
saveRDS(gsea_results,"/cluster/home/xyzhang_jh/projects/Xenograft/pig/GSEA_sscrofa.rds")

#node70
gsea_results_kegg <- readRDS(file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/pig/GSEA_sscrofa_kegg.rds")

pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/pig/facet_by_sign_kegg.pdf",width = 8,height = 10)
dotplot(gsea_results_kegg,split=".sign")+facet_grid(~.sign)+
  theme(
    axis.text = element_text(size = 8),  # 坐标轴文本大小
    axis.title = element_text(size = 12), # 坐标轴标题大小
    strip.text = element_text(size = 10), # 分面标签大小
    plot.title = element_text(size = 14)  # 图表标题大小
  ) +
  ggtitle("GSEA KEGG Analysis LT_vs_LC")  # 添加标题
dev.off()

gsea_results_kegg2 <- as.data.frame(gsea_results_kegg)
convert_genes <- function(gene_ids) {
  ids <- unlist(strsplit(gene_ids, "/"))
  gene_info <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Ss.eg.db")
}
# 应用转换函数到 core_enrichment 列
a <- sapply(gsea_results_kegg2$core_enrichment, convert_genes)

gsea_results_kegg2$core_enrichment <- sapply(a[2,], function(x) {
  paste(unlist(x), collapse = "/")
})
write.csv(gsea_results_kegg2, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/pig/gsea_kegg_res.csv")



# 运行 GSEA
gsea_results_go <- gseGO(geneList = gene_list,
                         OrgDb = "org.Ss.eg.db",  # 猪
                         nPerm = 1000,      # 随机化次数
                         minGSSize = 10,    # 基因集最小大小
                         maxGSSize = 500,   # 基因集最大大小
                         pvalueCutoff = 1,
                         verbose = TRUE)
saveRDS(gsea_results_go,"/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/pig/GSEA_pig_go.rds")
gsea_results_go <- readRDS("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/pig/GSEA_pig_go.rds")
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/pig/facet_by_sign_go.pdf",width = 8,height = 10)
dotplot(gsea_results_go,split=".sign")+facet_grid(~.sign)+
  theme(
    axis.text = element_text(size = 8),  # 坐标轴文本大小
    axis.title = element_text(size = 12), # 坐标轴标题大小
    strip.text = element_text(size = 10), # 分面标签大小
    plot.title = element_text(size = 14)  # 图表标题大小
  ) +
  ggtitle("GSEA GO Analysis LT_vs_LC")  # 添加标题
dev.off()

gsea_results_go2 <- as.data.frame(gsea_results_go)
convert_genes <- function(gene_ids) {
  ids <- unlist(strsplit(gene_ids, "/"))
  gene_info <- bitr(ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Ss.eg.db")
}
# 应用转换函数到 core_enrichment 列
a <- sapply(gsea_results_go2$core_enrichment, convert_genes)

gsea_results_go2$core_enrichment <- sapply(a[2,], function(x) {
  paste(unlist(x), collapse = "/")
})
write.csv(gsea_results_go2, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/GSEA/pig/gsea_go_res.csv")
















