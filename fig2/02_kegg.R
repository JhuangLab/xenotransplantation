.libPaths(new = c(.libPaths(), "/cluster/home/danyang_jh/sbin/R/R-4.3.0", "/cluster/home/danyang_jh/sbin/R/R-4.2.1","/cluster/home/xyzhang_jh/sbin/R"))
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "QFeatures", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "MSnbase", 
          "enrichplot", "org.Hs.eg.db", "parallel", "jhuanglabRNAseq", "jhuanglabGO", "limma", 
          "readxl", "writexl", "ComplexHeatmap", "circlize", "clusterProfiler", "MsCoreUtils", 
          "paletteer", "factoextra", "FactoMineR", "openxlsx","gdata","ggrepel","clusterProfiler","org.Ss.eg.db","DOSE","org.Mmu.eg.db")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project <- "liver"
dataset <- "wangrongrong"
species <- "Xeno_pig_monkey"
workdir <- glue::glue("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/KEGG")
workdir %>% fs::dir_create() %>% setwd()

# kegg monkey
fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
fc_table$gene_name <- rownames(fc_table)


gene_entrez <- bitr(fc_table$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mmu.eg.db")
df_Mmu_sig <- merge(fc_table, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL") %>%
  filter(FDR < 0.05 & abs(logFC) > 1)
saveRDS(df_Mmu_sig,"df_sig_mmu_FC1.rds")
df_Mmu_sig <- read_rds("df_sig_mmu_FC1.rds")
kegg_mmul <- enrichKEGG(gene = df_Mmu_sig$ENTREZID, 
                        organism = "mcc", pvalueCutoff = 0.05)
kegg_mmul_up <- enrichKEGG(gene = df_Mmu_sig[df_Mmu_sig$logFC > 0, ]$ENTREZID, 
                           organism = "mcc", pvalueCutoff = 0.05) 
kegg_mmul_down <- enrichKEGG(gene = df_Mmu_sig[df_Mmu_sig$logFC < 0, ]$ENTREZID, 
                             organism = "mcc", pvalueCutoff = 0.05)
gene_entrez <- df_Mmu_sig[,c("gene_name","ENTREZID")]

dt_kegg_mmul <- as.data.frame(kegg_mmul) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_kegg_mmul_up <- as.data.frame(kegg_mmul_up) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_kegg_mmul_down <- as.data.frame(kegg_mmul_down) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
fwrite(dt_kegg_mmul, "results/diff/FC_1/kegg_mmul_all.csv")
fwrite(dt_kegg_mmul_up, "results/diff/FC_1/kegg_mmul_up.csv")
fwrite(dt_kegg_mmul_down, "results/diff/FC_1/kegg_mmul_down.csv")

p1 <- dotplot(kegg_mmul) + ggtitle("KEGG all")
p2 <- dotplot(kegg_mmul_up) + ggtitle("KEGG up in post")
p3 <- dotplot(kegg_mmul_down) + ggtitle("KEGG down in post")
ggsave("results/diff/FC_1/dotplot_kegg_all_mmul.pdf", p1, width = 7, height = 8)
ggsave("results/diff/FC_1/dotplot_kegg_up_mmul.pdf", p2, width = 7, height = 6)
ggsave("results/diff/FC_1/dotplot_kegg_down_mmul.pdf", p3, width = 7, height = 7)

# kegg pig
fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
fc_table$gene_name <- rownames(fc_table)
gene_entrez <- bitr(fc_table$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Ss.eg.db")
df_ssc_sig <- merge(fc_table, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL") %>%
  filter(FDR < 0.05 & abs(logFC) > 1)
saveRDS(df_ssc_sig,"df_sig_ssc_FC1.rds")

df_ssc_sig <- read_rds("df_sig_ssc_FC1.rds")



kegg_ssc <- enrichKEGG(gene = df_ssc_sig$ENTREZID, 
                       organism = "ssc", pvalueCutoff = 0.05)
kegg_ssc_up <- enrichKEGG(gene = df_ssc_sig[df_ssc_sig$logFC > 0, ]$ENTREZID, 
                          organism = "ssc", pvalueCutoff = 0.05) 
kegg_ssc_down <- enrichKEGG(gene = df_ssc_sig[df_ssc_sig$logFC < 0, ]$ENTREZID, 
                            organism = "ssc", pvalueCutoff = 0.05)
gene_entrez <- df_ssc_sig[,c("gene_name","ENTREZID")]

dt_kegg_ssc <- as.data.frame(kegg_ssc) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_kegg_ssc_up <- as.data.frame(kegg_ssc_up) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_kegg_ssc_down <- as.data.frame(kegg_ssc_down) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
fwrite(dt_kegg_ssc, "results/diff/FC_1/kegg_ssc_all.csv")
fwrite(dt_kegg_ssc_up, "results/diff/FC_1/kegg_ssc_up.csv")
fwrite(dt_kegg_ssc_down, "results/diff/FC_1/kegg_ssc_down.csv")

p1 <- dotplot(kegg_ssc) + ggtitle("KEGG all")
p2 <- dotplot(kegg_ssc_up) + ggtitle("KEGG up in post")
p3 <- dotplot(kegg_ssc_down) + ggtitle("KEGG down in post")
ggsave("results/diff/FC_1/dotplot_kegg_all_ssc.pdf", p1, width = 7, height = 8)
ggsave("results/diff/FC_1/dotplot_kegg_up_ssc.pdf", p2, width = 7, height = 6)
ggsave("results/diff/FC_1/dotplot_kegg_down_ssc.pdf", p3, width = 7, height = 7)



##FC=0
# kegg monkey
fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
fc_table$gene_name <- rownames(fc_table)

gene_entrez <- bitr(fc_table$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mmu.eg.db")
df_Mmu_sig <- merge(fc_table, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL") %>%
  filter(FDR < 0.05)
saveRDS(df_Mmu_sig,"df_sig_mmu_FC0.rds")
df_Mmu_sig <- read_rds("df_sig_mmu_FC0.rds")
df_Mmu_sig$group <- ifelse(df_Mmu_sig$logFC > 0 ,"up","down") 

kegg_mmul <- enrichKEGG(gene = df_Mmu_sig$ENTREZID, 
                        organism = "mcc", pvalueCutoff = 0.05)
kegg_mmul_up <- enrichKEGG(gene = df_Mmu_sig[df_Mmu_sig$logFC > 0, ]$ENTREZID, 
                           organism = "mcc", pvalueCutoff = 0.05) 
kegg_mmul_down <- enrichKEGG(gene = df_Mmu_sig[df_Mmu_sig$logFC < 0, ]$ENTREZID, 
                             organism = "mcc", pvalueCutoff = 0.05)
gene_entrez <- df_Mmu_sig[,c("gene_name","ENTREZID")]

dt_kegg_mmul <- as.data.frame(kegg_mmul) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_kegg_mmul_up <- as.data.frame(kegg_mmul_up) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_kegg_mmul_down <- as.data.frame(kegg_mmul_down) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
fwrite(dt_kegg_mmul, "results/diff/FC_0/kegg_mmul_all.csv")
fwrite(dt_kegg_mmul_up, "results/diff/FC_0/kegg_mmul_up.csv")
fwrite(dt_kegg_mmul_down, "results/diff/FC_0/kegg_mmul_down.csv")

p1 <- dotplot(kegg_mmul) + ggtitle("KEGG all")
p2 <- dotplot(kegg_mmul_up) + ggtitle("KEGG up in post")
p3 <- dotplot(kegg_mmul_down) + ggtitle("KEGG down in post")
ggsave("results/diff/FC_0/dotplot_kegg_all_mmul.pdf", p1, width = 7, height = 8)
ggsave("results/diff/FC_0/dotplot_kegg_up_mmul.pdf", p2, width = 7, height = 6)
ggsave("results/diff/FC_0/dotplot_kegg_down_mmul.pdf", p3, width = 7, height = 7)

# kegg pig
fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
fc_table$gene_name <- rownames(fc_table)
gene_entrez <- bitr(fc_table$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Ss.eg.db")
df_ssc_sig <- merge(fc_table, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL") %>%
  filter(FDR < 0.05)
saveRDS(df_ssc_sig,"df_sig_ssc_FC0.rds")

df_ssc_sig <- read_rds("df_sig_ssc_FC0.rds")



kegg_ssc <- enrichKEGG(gene = df_ssc_sig$ENTREZID, 
                       organism = "ssc", pvalueCutoff = 0.05)
kegg_ssc_up <- enrichKEGG(gene = df_ssc_sig[df_ssc_sig$logFC > 0, ]$ENTREZID, 
                          organism = "ssc", pvalueCutoff = 0.05) 
kegg_ssc_down <- enrichKEGG(gene = df_ssc_sig[df_ssc_sig$logFC < 0, ]$ENTREZID, 
                            organism = "ssc", pvalueCutoff = 0.05)
gene_entrez <- df_ssc_sig[,c("gene_name","ENTREZID")]

dt_kegg_ssc <- as.data.frame(kegg_ssc) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_kegg_ssc_up <- as.data.frame(kegg_ssc_up) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_kegg_ssc_down <- as.data.frame(kegg_ssc_down) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
fwrite(dt_kegg_ssc, "results/diff/FC_0/kegg_ssc_all.csv")
fwrite(dt_kegg_ssc_up, "results/diff/FC_0/kegg_ssc_up.csv")
fwrite(dt_kegg_ssc_down, "results/diff/FC_0/kegg_ssc_down.csv")

p1 <- dotplot(kegg_ssc) + ggtitle("KEGG all")
p2 <- dotplot(kegg_ssc_up) + ggtitle("KEGG up in post")
p3 <- dotplot(kegg_ssc_down) + ggtitle("KEGG down in post")
ggsave("results/diff/FC_0/dotplot_kegg_all_ssc.pdf", p1, width = 7, height = 8)
ggsave("results/diff/FC_0/dotplot_kegg_up_ssc.pdf", p2, width = 7, height = 6)
ggsave("results/diff/FC_0/dotplot_kegg_down_ssc.pdf", p3, width = 7, height = 7)


##绘制指定的通路
##ssc up
selected_pathways_ssc <- c("Phagosome","Efferocytosis","ECM-receptor interaction","Complement and coagulation cascades","Glycolysis / Gluconeogenesis","Biosynthesis of amino acids","Spliceosome","Thermogenesis","Bile secretion","Non-alcoholic fatty liver disease")
##mmul uo
selected_pathways_mmul <- c("Citrate cycle (TCA cycle)","Fatty acid degradation","Oxidative phosphorylation","Phagosome","Efferocytosis","Focal adhesion","ECM-receptor interaction","Complement and coagulation cascades","Platelet activation","Non-alcoholic fatty liver disease")




dt_kegg_ssc_up <- read.csv(file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/KEGG/results/diff/FC_0/kegg_ssc_up.csv")
dt_kegg_mmul_up <- read.csv(file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/KEGG/results/diff/FC_0/kegg_mmul_up.csv")
p2 <- dotplot(dt_kegg_ssc_up) + ggtitle("KEGG up in post")
p3 <- dotplot(dt_kegg_mmul_up) + ggtitle("KEGG down in post")

df_Mmu_sig <- read_rds("df_sig_mmu_FC0.rds")
kegg_mmul_up <- enrichKEGG(gene = df_Mmu_sig[df_Mmu_sig$logFC > 0, ]$ENTREZID, 
                           organism = "mcc", pvalueCutoff = 0.05)

kk_filtered <- kegg_mmul_up@result %>%
  filter(Description %in% selected_pathways_mmul)
kk_filtered <- new("enrichResult", result = kk_filtered)
p1 <- dotplot(kk_filtered, showCategory = length(selected_pathways_mmul))
ggsave(plot = p1, "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/KEGG/results/diff/FC_0/dotplot_kegg_up_mmul_selected.pdf", width = 7, height = 6)



df_ssc_sig <- read_rds("df_sig_ssc_FC0.rds")
kegg_ssc_up <- enrichKEGG(gene = df_ssc_sig[df_ssc_sig$logFC > 0, ]$ENTREZID, 
                          organism = "ssc", pvalueCutoff = 0.05) 
kk_filtered <- kegg_ssc_up@result %>%
  filter(Description %in% selected_pathways_ssc)
kk_filtered <- new("enrichResult", result = kk_filtered)
p1 <- dotplot(kk_filtered, showCategory = length(selected_pathways_ssc))
ggsave(plot = p1, "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/KEGG/results/diff/FC_0/dotplot_kegg_up_ssc_selected.pdf", width = 7, height = 6)













