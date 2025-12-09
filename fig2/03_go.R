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
workdir <- glue::glue("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/GO")
workdir %>% fs::dir_create() %>% setwd()


fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
fc_table$gene_name <- rownames(fc_table)

# go monkey
gene_entrez <- bitr(fc_table$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mmu.eg.db")
df_Mmu_sig <- merge(fc_table, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL") %>%
  filter(FDR < 0.05 & abs(logFC) > 1)
saveRDS(df_Mmu_sig,"df_sig_mmu.rds")
df_Mmu_sig <- read_rds("df_sig_mmu.rds")
go_mmul <- enrichGO(gene = df_Mmu_sig$ENTREZID, 
                    OrgDb = "org.Mmu.eg.db", 
                    pvalueCutoff = 0.05)
go_mmul_up <- enrichGO(gene = df_Mmu_sig[df_Mmu_sig$logFC > 0, ]$ENTREZID, 
                           OrgDb = "org.Mmu.eg.db",  pvalueCutoff = 0.05) 
go_mmul_down <- enrichGO(gene = df_Mmu_sig[df_Mmu_sig$logFC < 0, ]$ENTREZID, 
                             OrgDb = "org.Mmu.eg.db",  pvalueCutoff = 0.05)
gene_entrez <- df_Mmu_sig[,c("gene_name","ENTREZID")]

dt_go_mmul <- as.data.frame(go_mmul) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_go_mmul_up <- as.data.frame(go_mmul_up) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_go_mmul_down <- as.data.frame(go_mmul_down) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
fwrite(dt_go_mmul, "results/diff/FC_1/go_mmul_all.csv")
fwrite(dt_go_mmul_up, "results/diff/FC_1/go_mmul_up.csv")
fwrite(dt_go_mmul_down, "results/diff/FC_1/go_mmul_down.csv")

p1 <- dotplot(go_mmul) + ggtitle("go all")
p2 <- dotplot(go_mmul_up) + ggtitle("go up in post")
p3 <- dotplot(go_mmul_down) + ggtitle("go down in post")
ggsave("results/diff/FC_1/dotplot_go_all_mmul.pdf", p1, width = 7, height = 8)
ggsave("results/diff/FC_1/dotplot_go_up_mmul.pdf", p2, width = 7, height = 6)
ggsave("results/diff/FC_1/dotplot_go_down_mmul.pdf", p3, width = 7, height = 7)

# go pig
fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
fc_table$gene_name <- rownames(fc_table)
gene_entrez <- bitr(fc_table$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Ss.eg.db")
df_ssc_sig <- merge(fc_table, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL") %>%
  filter(FDR < 0.05 & abs(logFC) > 1)
saveRDS(df_ssc_sig,"df_sig_ssc.rds")

df_ssc_sig <- read_rds("df_sig_ssc.rds")



go_ssc <- enrichGO(gene = df_ssc_sig$ENTREZID, 
                   OrgDb = "org.Ss.eg.db", pvalueCutoff = 0.05)
go_ssc_up <- enrichGO(gene = df_ssc_sig[df_ssc_sig$logFC > 0, ]$ENTREZID, 
                      OrgDb = "org.Ss.eg.db", pvalueCutoff = 0.05) 
go_ssc_down <- enrichGO(gene = df_ssc_sig[df_ssc_sig$logFC < 0, ]$ENTREZID, 
                        OrgDb = "org.Ss.eg.db", pvalueCutoff = 0.05)
gene_entrez <- df_ssc_sig[,c("gene_name","ENTREZID")]

dt_go_ssc <- as.data.frame(go_ssc) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_go_ssc_up <- as.data.frame(go_ssc_up) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_go_ssc_down <- as.data.frame(go_ssc_down) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
fwrite(dt_go_ssc, "results/diff/FC_1/go_ssc_all.csv")
fwrite(dt_go_ssc_up, "results/diff/FC_1/go_ssc_up.csv")
fwrite(dt_go_ssc_down, "results/diff/FC_1/go_ssc_down.csv")

p1 <- dotplot(go_ssc) + ggtitle("go all")
p2 <- dotplot(go_ssc_up) + ggtitle("go up in post")
p3 <- dotplot(go_ssc_down) + ggtitle("go down in post")
ggsave("results/diff/FC_1/dotplot_go_all_ssc.pdf", p1, width = 7, height = 8)
ggsave("results/diff/FC_1/dotplot_go_up_ssc.pdf", p2, width = 7, height = 6)
ggsave("results/diff/FC_1/dotplot_go_down_ssc.pdf", p3, width = 7, height = 7)




## FC = 0
fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
fc_table$gene_name <- rownames(fc_table)

# go monkey
gene_entrez <- bitr(fc_table$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mmu.eg.db")
df_Mmu_sig <- merge(fc_table, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL") %>%
  filter(FDR < 0.05)
saveRDS(df_Mmu_sig,"df_sig_mmu.rds")
df_Mmu_sig <- read_rds("df_sig_mmu.rds")
go_mmul <- enrichGO(gene = df_Mmu_sig$ENTREZID, 
                    OrgDb = "org.Mmu.eg.db", 
                    pvalueCutoff = 0.05)
go_mmul_up <- enrichGO(gene = df_Mmu_sig[df_Mmu_sig$logFC > 0, ]$ENTREZID, 
                       OrgDb = "org.Mmu.eg.db",  pvalueCutoff = 0.05) 
go_mmul_down <- enrichGO(gene = df_Mmu_sig[df_Mmu_sig$logFC < 0, ]$ENTREZID, 
                         OrgDb = "org.Mmu.eg.db",  pvalueCutoff = 0.05)
gene_entrez <- df_Mmu_sig[,c("gene_name","ENTREZID")]

dt_go_mmul <- as.data.frame(go_mmul) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_go_mmul_up <- as.data.frame(go_mmul_up) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_go_mmul_down <- as.data.frame(go_mmul_down) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
fwrite(dt_go_mmul, "results/diff/FC_0/go_mmul_all.csv")
fwrite(dt_go_mmul_up, "results/diff/FC_0/go_mmul_up.csv")
fwrite(dt_go_mmul_down, "results/diff/FC_0/go_mmul_down.csv")

p1 <- dotplot(go_mmul) + ggtitle("go all")
p2 <- dotplot(go_mmul_up) + ggtitle("go up in post")
p3 <- dotplot(go_mmul_down) + ggtitle("go down in post")
ggsave("results/diff/FC_0/dotplot_go_all_mmul.pdf", p1, width = 7, height = 8)
ggsave("results/diff/FC_0/dotplot_go_up_mmul.pdf", p2, width = 7, height = 6)
ggsave("results/diff/FC_0/dotplot_go_down_mmul.pdf", p3, width = 7, height = 7)

# go pig
fc_table <- read.csv("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/single-sample-deg-result_0.1.csv", row.names = 1)
fc_table$gene_name <- rownames(fc_table)
gene_entrez <- bitr(fc_table$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Ss.eg.db")
df_ssc_sig <- merge(fc_table, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL") %>%
  filter(FDR < 0.05)
saveRDS(df_ssc_sig,"df_sig_ssc.rds")

df_ssc_sig <- read_rds("df_sig_ssc.rds")



go_ssc <- enrichGO(gene = df_ssc_sig$ENTREZID, 
                   OrgDb = "org.Ss.eg.db", pvalueCutoff = 0.05)
go_ssc_up <- enrichGO(gene = df_ssc_sig[df_ssc_sig$logFC > 0, ]$ENTREZID, 
                      OrgDb = "org.Ss.eg.db", pvalueCutoff = 0.05) 
go_ssc_down <- enrichGO(gene = df_ssc_sig[df_ssc_sig$logFC < 0, ]$ENTREZID, 
                        OrgDb = "org.Ss.eg.db", pvalueCutoff = 0.05)
gene_entrez <- df_ssc_sig[,c("gene_name","ENTREZID")]

dt_go_ssc <- as.data.frame(go_ssc) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_go_ssc_up <- as.data.frame(go_ssc_up) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
dt_go_ssc_down <- as.data.frame(go_ssc_down) %>% 
  tidyr::separate_rows(geneID, sep = '/') %>% as.data.table %>%
  merge(gene_entrez, all.x = TRUE, by.x = "geneID", by.y = "ENTREZID") %>%
  .[order(ID), ]
fwrite(dt_go_ssc, "results/diff/FC_0/go_ssc_all.csv")
fwrite(dt_go_ssc_up, "results/diff/FC_0/go_ssc_up.csv")
fwrite(dt_go_ssc_down, "results/diff/FC_0/go_ssc_down.csv")

p1 <- dotplot(go_ssc) + ggtitle("go all")
p2 <- dotplot(go_ssc_up) + ggtitle("go up in post")
p3 <- dotplot(go_ssc_down) + ggtitle("go down in post")
ggsave("results/diff/FC_0/dotplot_go_all_ssc.pdf", p1, width = 7, height = 8)
ggsave("results/diff/FC_0/dotplot_go_up_ssc.pdf", p2, width = 7, height = 6)
ggsave("results/diff/FC_0/dotplot_go_down_ssc.pdf", p3, width = 7, height = 7)