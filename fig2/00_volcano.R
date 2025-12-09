.libPaths(new = c(.libPaths(), "/cluster/home/danyang_jh/sbin/R/R-4.3.0", "/cluster/home/danyang_jh/sbin/R/R-4.2.1","/cluster/home/xyzhang_jh/sbin/R"))
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "QFeatures", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "MSnbase", "msqrob2", 
          "enrichplot", "org.Hs.eg.db", "parallel", "jhuanglabRNAseq", "jhuanglabGO", "limma", 
          "readxl", "writexl", "ComplexHeatmap", "circlize", "clusterProfiler", "MsCoreUtils", 
          "missForest", "paletteer", "factoextra", "FactoMineR", "openxlsx","gdata","ggrepel","edgeR")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "liver"
dataset <- "wangrongrong"
species <- "Xeno_pig_monkey"
workdir <- glue::glue("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein")
workdir %>% fs::dir_create() %>% setwd()

install.packages("/cluster/home/xyzhang_jh/sbin/R/org.Mmu.eg.db_3.21.0.tar.gz", repos = NULL, type = "source", lib = "/cluster/home/xyzhang_jh/sbin/R")

my_theme1 <- theme_classic(base_size = 8) + 
  theme(legend.position = "right", legend.key.size = unit(3, "mm"), 
        plot.title = element_text(hjust = .5), axis.line = element_blank(), 
        panel.border = element_rect(fill = NA, linewidth = .5))

output_dir <- "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein" %>% checkdir()


all_exp <- read.delim("/cluster/home/jhuang/projects/liver/data/zhangwei/pig/protein/report/report-liver/data/3.comparison/PIT_CON-IT/PIT_CON-IT.all.xls")
all_exp <- all_exp[,c("gene_name","PIT","CON.IT")]
all_exp <- all_exp[!grepl("-", all_exp$gene_name), ] %>% na.omit()
all_exp <- all_exp[!duplicated(all_exp$gene_name),]
rownames(all_exp) <- NULL

all_exp <- all_exp %>% column_to_rownames(.,"gene_name")

# 应用 log2(exp + 1) 转换
#all_exp <- log2(all_exp + 1)
#define the sample groups
group <- c("post", "control")
#edgeR
sample_info <- data_frame(
  Sample = colnames(all_exp), 
  Group = group
)
all_exp <- as.matrix(all_exp)
dge <- DGEList(counts = all_exp, group = group)
dge <- calcNormFactors(dge)

bcv = 0.1#设置bcv为0.1
et <- exactTest(dge, dispersion=bcv^2)
DEG_edgeR <- as.data.frame(topTags(et, n = nrow(dge$counts)))

DEG_edgeR[DEG_edgeR$FDR < 0.05, ] %>% dim()


write.csv(DEG_edgeR, '/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv', quote = FALSE)





dge <- DGEList(counts = all_exp, group = group)
dge <- calcNormFactors(dge)

bcv = 0.4#设置bcv为0.4
et <- exactTest(dge, dispersion=bcv^2)
DEG_edgeR <- as.data.frame(topTags(et, n = nrow(dge$counts)))
head(DEG_edgeR)
DEG_edgeR[DEG_edgeR$FDR < 0.05, ] %>% dim()

write.csv(DEG_edgeR, '/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.4.csv', quote = FALSE)


#volcanomap
results <- read.csv('/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv', row.names = 1)
output_dir <- "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein"
comp <- "post_vs_control"
results <- results %>% as.data.frame() %>% 
  mutate(trend = case_when((FDR < 0.05) & (logFC >= 0.5) ~ "up", 
                           (FDR < 0.05) & (logFC <= -0.5) ~ "down", 
                           TRUE ~ "Not-sig"))
results$gene_name <- rownames(results)
results <- results[!duplicated(results$gene_name),]
top_diff_genes <- results[order(results$FDR), ][1:10, ] %>% .[.$trend != "Not-sig", "gene_name"]

label_frame <- data_frame(gene_name = c(unique(top_diff_genes)),
                          tag = c(unique(top_diff_genes)))
results <- left_join(results, label_frame, by = "gene_name")
write.csv(results, file = glue(output_dir,"/single_volcano_filtered_padj_mmul.csv"))



vol1 <- ggplot(results, aes(x = logFC, y = -log10(FDR), color = trend, fill = trend)) +
  ylim(0, max(-log10(results$FDR))+0.5)+
  #xlim(-2,2)+
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
ggsave(glue(output_dir,"/volcano_padj_LT_vs_LC.pdf"), vol1, width = 5, height = 5)







