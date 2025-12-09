suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(rlang)
  library(ggplot2)
  library(edgeR)
  library(tximport)
  library(clusterProfiler)
  library(data.table)
})

wdir <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/"
setwd(wdir)

# import salmon quant files
samples <- data.frame(condition = factor(c("LC", "LT")),
                      row.names = c("LC", "LT"))
files <- file.path(wdir, "results/salmon_quant", c("LC", "LT"), "quant.sf")
tx2gene <- read.delim("ref/combined/tx2gene_combined.tsv", header = FALSE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# save raw counts & fpkm
mat_ct <- txi$counts
colnames(mat_ct) <- c("LC", "LT")
write.table(mat_ct, "results/salmon_quant/counts.txt")
mat_fpkm <- rpkm(mat_ct, gene.length = txi$length)
write.table(mat_fpkm, "results/salmon_quant/counts_fpkm.txt")

mat_ct_mk <- mat_ct[grep("Mmul_", rownames(mat_ct)), ]
rownames(mat_ct_mk) <- sub("Mmul_", "", rownames(mat_ct_mk))
write.table(mat_ct_mk, "results/salmon_quant/counts_mmul.txt")
mat_fpkm_mk <- mat_fpkm[grep("Mmul_", rownames(mat_ct)), ]
rownames(mat_fpkm_mk) <- sub("Mmul_", "", rownames(mat_fpkm_mk))
write.table(mat_fpkm_mk, "results/salmon_quant/counts_fpkm_mmul.txt")

mat_ct_pg <- mat_ct[grep("Sscrofa_", rownames(mat_ct)), ]
rownames(mat_ct_pg) <- sub("Sscrofa_", "", rownames(mat_ct_pg))
write.table(mat_ct_pg, "results/salmon_quant/counts_sscrofa.txt")
mat_fpkm_pg <- mat_fpkm[grep("Sscrofa_", rownames(mat_ct)), ]
rownames(mat_fpkm_pg) <- sub("Sscrofa_", "", rownames(mat_fpkm_pg))
write.table(mat_fpkm_pg, "results/salmon_quant/counts_fpkm_sscrofa.txt")

# edgeR
group <- samples$condition
y <- DGEList(counts = txi$counts, group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y)
et <- exactTest(y, dispersion = 0.1) # 0.1 for data on genetically identical model organisms when no replicates

# diff genes
df_diff <- et$table
df_diff$fdr <- p.adjust(df_diff$PValue, method = "fdr")
write.csv(df_diff, "results/diff/diff_genes_combined.csv")
df_sig <- df_diff %>%
  filter(fdr < 0.05 & abs(logFC) > 1) %>%
  arrange(fdr) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  left_join(tibble::rownames_to_column(as.data.frame(mat_fpkm), var = "gene_name"))
write.csv(df_sig, "results/diff/diff_genes_combined_sig.csv", row.names = FALSE)

# separate monkey pig
df_mmul <- tibble::rownames_to_column(df_diff, var = "gene_name") %>%
  filter(stringr::str_detect(gene_name, "Mmul_")) %>%
  mutate(gene_name = sub("Mmul_", "", gene_name))
df_sscrofa <- tibble::rownames_to_column(df_diff, var = "gene_name") %>%
  filter(stringr::str_detect(gene_name, "Sscrofa_")) %>%
  mutate(gene_name = sub("Sscrofa_", "", gene_name))
write.csv(df_mmul, "results/diff/diff_genes_mmul.csv", row.names = FALSE)
write.csv(df_sscrofa, "results/diff/diff_genes_sscrofa.csv", row.names = FALSE)

# kegg monkey
df_mmul <- read.csv("results/diff/diff_genes_mmul.csv")
gene_entrez <- bitr(df_mmul$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mmu.eg.db")
df_mmul_sig <- merge(df_mmul, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL") %>%
  filter(fdr < 0.05 & abs(logFC) > 1)
kegg_mmul <- enrichKEGG(gene = df_mmul_sig$ENTREZID, 
                   organism = "mcc", pvalueCutoff = 0.05)
kegg_mmul_up <- enrichKEGG(gene = df_mmul_sig[df_mmul_sig$logFC > 0, ]$ENTREZID, 
                      organism = "mcc", pvalueCutoff = 0.05) 
kegg_mmul_down <- enrichKEGG(gene = df_mmul_sig[df_mmul_sig$logFC < 0, ]$ENTREZID, 
                        organism = "mcc", pvalueCutoff = 0.05)

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
fwrite(dt_kegg_mmul, "results/diff/kegg_mmul_all.csv")
fwrite(dt_kegg_mmul_up, "results/diff/kegg_mmul_up.csv")
fwrite(dt_kegg_mmul_down, "results/diff/kegg_mmul_down.csv")

p1 <- dotplot(kegg_mmul) + ggtitle("KEGG all")
p2 <- dotplot(kegg_mmul_up) + ggtitle("KEGG up in LT")
p3 <- dotplot(kegg_mmul_down) + ggtitle("KEGG down in LT")
ggsave("results/diff/dotplot_kegg_all_mmul.pdf", p1, width = 7, height = 8)
ggsave("results/diff/dotplot_kegg_up_mmul.pdf", p2, width = 7, height = 6)
ggsave("results/diff/dotplot_kegg_down_mmul.pdf", p3, width = 7, height = 7)

# kegg pig
df_ssc <- read.csv("results/diff/diff_genes_sscrofa.csv")
gene_entrez <- bitr(df_ssc$gene_name, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Ss.eg.db")
df_ssc_sig <- merge(df_ssc, gene_entrez, all.x = TRUE, by.x = "gene_name", by.y = "SYMBOL") %>%
  filter(fdr < 0.05 & abs(logFC) > 1)
kegg_ssc <- enrichKEGG(gene = df_ssc_sig$ENTREZID, 
                        organism = "ssc", pvalueCutoff = 0.05)
kegg_ssc_up <- enrichKEGG(gene = df_ssc_sig[df_ssc_sig$logFC > 0, ]$ENTREZID, 
                           organism = "ssc", pvalueCutoff = 0.05) 
kegg_ssc_down <- enrichKEGG(gene = df_ssc_sig[df_ssc_sig$logFC < 0, ]$ENTREZID, 
                             organism = "ssc", pvalueCutoff = 0.05)

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
fwrite(dt_kegg_ssc, "results/diff/kegg_ssc_all.csv")
fwrite(dt_kegg_ssc_up, "results/diff/kegg_ssc_up.csv")
fwrite(dt_kegg_ssc_down, "results/diff/kegg_ssc_down.csv")

p1 <- dotplot(kegg_ssc) + ggtitle("KEGG all")
p2 <- dotplot(kegg_ssc_up) + ggtitle("KEGG up in LT")
p3 <- dotplot(kegg_ssc_down) + ggtitle("KEGG down in LT")
ggsave("results/diff/dotplot_kegg_all_ssc.pdf", p1, width = 7, height = 8)
ggsave("results/diff/dotplot_kegg_up_ssc.pdf", p2, width = 7, height = 6)
ggsave("results/diff/dotplot_kegg_down_ssc.pdf", p3, width = 7, height = 7)
