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

wdir <- "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/edgeR_test"
setwd(wdir)

# import salmon quant files
all_exp <- read.delim("/cluster/home/jhuang/projects/liver/data/zhangwei/pig/protein/report/report-liver/data/3.comparison/PIT_CON-IT/PIT_CON-IT.all.xls")
all_exp <- all_exp[,c("gene_name","PIT","CON.IT")]
all_exp <- all_exp[!grepl("-", all_exp$gene_name), ] %>% na.omit()
all_exp <- all_exp[!duplicated(all_exp$gene_name),]
rownames(all_exp) <- NULL

all_exp <- all_exp %>% column_to_rownames(.,"gene_name")

# save raw counts & fpkm
samples <- data.frame(condition = factor(c("LT", "LC")),
                      row.names = c("LT", "LC"))


# edgeR
group <- samples$condition
y <- DGEList(counts = all_exp, group = group)
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- normLibSizes(y)
et <- exactTest(y, dispersion = 0.1) # 0.1 for data on genetically identical model organisms when no replicates

# diff genes
df_diff <- et$table
df_diff$fdr <- p.adjust(df_diff$PValue, method = "fdr")
write.csv(df_diff, "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/edgeR_test/diff_genes_combined.csv")
df_sig <- df_diff %>%
  filter(fdr < 0.05 & abs(logFC) > 1) %>%
  arrange(fdr) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  left_join(tibble::rownames_to_column(as.data.frame(mat_fpkm), var = "gene_name"))
write.csv(df_sig, "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/protein/edgeR_test/results/diff/diff_genes_combined_sig.csv", row.names = FALSE)