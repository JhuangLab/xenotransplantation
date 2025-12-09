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
workdir <- glue::glue("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/WGCNA")
workdir %>% fs::dir_create() %>% setwd()

output_dir <- "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/WGCNA"


library(FSA, lib.loc = "/cluster/home/zhangjie_jh/my_Rlib")
library(dunn.test, lib.loc = "/cluster/home/zhangjie_jh/my_Rlib")



my_theme1 <- theme_classic(base_size = 8) + 
  theme(legend.position = "right", legend.key.size = unit(3, "mm"), 
        plot.title = element_text(hjust = .5), axis.line = element_blank(), 
        panel.border = element_rect(fill = NA, linewidth = .5))


file_path <- "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm.txt"
total_exp <- read.table(file_path)
names(total_exp) <- c("gene_name","logFC","logCPM","PValue","fdr")
# 去掉全是 0 的行
total_exp <- total_exp[rowSums(total_exp != 0) > 0, ]

## the expression data requires format -- rownames = samplename colnames = gene_name
#total_exp
log_frame <- data.frame(
  LC = log10(total_exp$LC + 1),
  LT = log10(total_exp$LT + 1)
)
rownames(log_frame) <- rownames(total_exp)

var_phos <- apply(total_exp, 1, var)
phosps <- var_phos > var_phos[names(sort(-var_phos)[2000])]
rnaData <- total_exp[phosps, ]


write.csv(rnaData, glue("{workdir}/clean_rnaseq_FPKM_code_gene_highvar2000_total.csv"))  
rnaData_t <- rnaData %>% t()

# don't log
# rnaData_log <- log2(rnaData_t + 1)
rnaData_log <- rnaData_t


#gsg = goodSamplesGenes(rnaData_log, verbose = 3)



# 异常值检测
sampleTree = hclust(dist(rnaData_log), method = "average")

## 异常tree检测
pdf(glue("{workdir}/tree.pdf"), width = 20, height = 9)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

pdf(glue("{workdir}/tree_check.pdf"), width = 20, height = 9)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
abline(h=235,col="red")
dev.off()
# #人工判断是否有离群值，进行剪枝操作
# clust <- cutreeStatic(sampleTree, cutHeight= 235)
# 
# keepSamples <- (clust==1)
# 用剪枝后数据进行后续分析
data = rnaData_log
write.csv(as.data.frame(data), glue("{workdir}/WGCNA_input_data.csv"))
# 自动选取软阈值，"power = sft$powerEstimate"，若“powerEstimate”为NA，需要手动确定阈值
# 手动选取阈值时R^2尽可能大，"Mean Connectivity"尽可能小

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(data, powerVector=powers,
                        networkType="unsigned", verbose=5)
saveRDS(sft, glue("{workdir}/WGCNA_select_softThread.rds"))

pdf(glue("{workdir}/select_softThread.pdf"), width = 12, height = 8)
par(mfrow = c(1,2))
# 理论此处应该大于等于0.85！！！！！
cex1 = 0.85
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),font.axis=3,cex.axis=2,font.lab=3,cex.lab=2)
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.85,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"),font.axis=3,cex.axis=2,font.lab=3,cex.lab=2)
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
dev.off()


power = sft$powerEstimate

net = blockwiseModules(data, power = power, maxBlockSize = ncol(data),
                       TOMType = "signed", minModuleSize = 20, randomSeed = 54321, 
                       numericLabels = TRUE, # reassignThreshold = 0,
                       pamRespectsDendro = FALSE, saveTOMs = TRUE, nThreads = 64, 
                       saveTOMFileBase = glue("{workdir}/blockwiseTOM"), verbose = 3)

saveRDS(net, glue("{workdir}/WGCNA_net.rds"))

# net <- readRDS(glue("{workdir}/Fig1/WGCNA/WGCNA_net.rds"))
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs_tmp = moduleEigengenes(data, moduleColors)
MEs = MEs_tmp$eigengenes
MEs = orderMEs(MEs)

gene_color_map <- data.frame(gene_name = colnames(data),
                             module_color = glue("ME{MEs_tmp$validColors}"))

write.csv(gene_color_map, glue("{workdir}/gene_color_map.csv"))

write.csv(MEs, glue("{workdir}/module_score.csv"))

color = moduleColors[net$blockGenes[[1]]]

pdf(glue("{workdir}/WGCNA_module.pdf"), width = 8, height = 6)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Cluster Dendrogram",
                    font.main=3,cex.main=2.5,
                    cex.colorLabels=1, font.axis=3,cex.axis=2,font.lab=3,cex.lab=2)
dev.off()


# barplot TNM grade-----------------------------------------------------------------
# 定义起始和结束颜色
color.bin <- c("#00599F", "#d80700")

# 创建渐变调色板函数
pal <- colorRampPalette(color.bin)

# 生成 3 个渐变颜色（包括两端）
group_col <- jhtools::show_me_the_colors(config_fn, "groups")

MEs$Patient_ID <- rownames(MEs)
names(sampinfo) <- c("Patient_ID","groups")
module_score <-  left_join(MEs,sampinfo, by = "Patient_ID") %>% 
  dplyr::select(Patient_ID, groups, everything())

plot.data <- module_score %>% 
  dplyr::select(-groups) %>% 
  column_to_rownames(var = "Patient_ID")

plot.info <- NULL
dunn_test_df <- NULL
module.name <- colnames(plot.data)
plot.stat <- data.frame(Module = module.name,
                        P = NA,
                        Padj = NA)

for (i in 1:ncol(plot.data)) {
  sub <- data.frame(Sample = rownames(plot.data),
                    Score = as.numeric(plot.data[, i]), 
                    ScoreScale = scale(as.numeric(plot.data[, i])),
                    Module = colnames(plot.data)[i],
                    Type = module_score$groups)
  
  # 使用 Kruskal-Wallis 检验（三分类变量）
  plot.stat$P[i] <- kruskal.test(ScoreScale ~ Type, data = sub)$p.value
  dunn_test_tmp <- dunnTest(ScoreScale ~ Type, data = sub, method = "bh")$res
  dunn_test_tmp$Module_name <- colnames(plot.data)[i]
  dunn_test_df <- rbind(dunn_test_df, dunn_test_tmp)
  
  plot.info <- rbind(plot.info, sub)
}
plot.info <- plot.info %>% mutate(Type=factor(as.character(Type), levels = c("Control", "Rejection-pre", "Rejection-post","ABOi"))) 
plot.stat <- plot.stat %>% mutate(Padj=p.adjust(P, method = "fdr"))


pdf(glue("{workdir}/module_score_in_groups.pdf"), width = 25)
p <- ggbarplot(
  plot.info,
  x = "Module",
  y = "ScoreScale",
  color = "Type",
  fill = "Type",
  palette = group_col,
  width = 0.5,
  size = 0,
  add = "mean_se",
  add.params = list(width = 0.5),
  order = module.name,
  position = position_dodge(0.6),
  xlab = "", 
  ylab = "Module Score of RNAseq"
) +
  theme_base() +
  geom_hline(yintercept = 0, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(plot.background = element_blank())

# 添加 Kruskal-Wallis 整体检验结果
p <- p + stat_compare_means(
  aes(group = Type), 
  method = "kruskal.test",
  label = "p.signif",
  label.y = 1.5,
  size = 3
)
print(p)
dev.off()

write.csv(plot.stat, glue("{workdir}/module_kruskal_test_results_groups.csv"))
write.csv(dunn_test_df, glue("{workdir}/module_dunn_test_TNM_groups.csv"))


# daizuo barplot survival grade---------------------------------------------------------

survival_state <- read_excel("/cluster/home/xyzhang_jh/projects/panlab/data/Survival_format_matrix.xlsx") %>% 
  mutate(Patient_ID = glue("P{META}"))
survival_state$surstat[survival_state$surstat == "1"] <- "alive"
survival_state$surstat[survival_state$surstat == "0"] <- "death"
survival_state$surstat <- factor(module_score$surstat, levels = c("death","alive"))

MEs <- read.csv("/cluster/home/xyzhang_jh/projects/panlab/analysis/baixueli/human/metabolome/wgcna_lc/module_score.csv", row.names = 1)
MEs$Patient_ID <- rownames(MEs)
MEs <- MEs %>% dplyr::filter(!rownames(MEs) %in% "P70") ## P70 with no sampleinfo
module_score <-  left_join(MEs,
                           survival_state[,c("surstat", "Patient_ID")], by = "Patient_ID") %>% 
  dplyr::select(Patient_ID, surstat, everything())
plot.data <- module_score %>% 
  dplyr::select(-surstat) %>% 
  column_to_rownames(var = "Patient_ID")

plot.info <- NULL
dunn_test_df <- NULL
module.name <- colnames(plot.data)
plot.stat <- data.frame(Module = module.name,
                        P = NA,
                        Padj = NA)

for (i in 1:ncol(plot.data)) {
  sub <- data.frame(Sample = rownames(plot.data),
                    Score = as.numeric(plot.data[, i]), 
                    ScoreScale = scale(as.numeric(plot.data[, i])),
                    Module = colnames(plot.data)[i],
                    Type = module_score$surstat)
  
  # 使用 Kruskal-Wallis 检验（三分类变量）
  plot.stat$P[i] <- kruskal.test(ScoreScale ~ Type, data = sub)$p.value
  
  dunn_test_tmp <- wilcox.test(ScoreScale ~ Type, data = sub, method = "bh")$res
  
  dunn_test_tmp$Module_name <- colnames(plot.data)[i]
  dunn_test_df <- rbind(dunn_test_df, dunn_test_tmp)
  plot.info <- rbind(plot.info, sub)
}


plot.info <- plot.info %>% mutate(Type=factor(as.character(Type), levels = c("0","1"))) 
plot.stat <- plot.stat %>% mutate(Padj=p.adjust(P, method = "fdr"))

pdf(glue("{workdir}/module_score_in_pathology_groups.pdf"),
    width = 25)
p <- ggbarplot(
  plot.info,
  x = "Module",
  y = "ScoreScale",
  color = "Type",
  fill = "Type",
  palette = colors_3,
  width = 0.5,
  size = 0,
  add = "mean_se",
  add.params = list(width = 0.5),
  order = module.name,
  position = position_dodge(0.6),
  xlab = "", 
  ylab = "Module Score of RNAseq"
) +
  theme_base() +
  geom_hline(yintercept = 0, color = "black") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(plot.background = element_blank())

# 添加 Kruskal-Wallis 整体检验结果
p <- p + stat_compare_means(
  aes(group = Type), 
  method = "kruskal.test",
  label = "p.signif",
  label.y = 1.5,
  size = 3
)
print(p)
dev.off()

write.csv(plot.stat, glue("{workdir}/module_kruskal_test_results_pathology.csv"))
write.csv(dunn_test_df, glue("{workdir}/module_dunn_test_pathology.csv"))



# 获取模块特征向量, 模块与表型相关联------------------------------------------------
MEs <- read.csv("module_score.csv", row.names = 1)

MEs_df <- MEs %>% 
  left_join(., sampinfo, by = "Patient_ID")

MEs_df$groups <- factor(MEs_df$groups, levels = c("Control", "Rejection-pre", "Rejection-post","ABOi"))

MEs_df <- MEs_df %>% column_to_rownames(var = "Patient_ID")


# 进行 One-Hot 编码
trait_encoded <- model.matrix(~ 0 + MEs_df$groups)  # 不加截距项
colnames(trait_encoded) <- c("Control", "Rejection-pre", "Rejection-post","ABOi")

MEs <- MEs %>% dplyr::select(-(Patient_ID))
modTraitCor <- WGCNA::cor(MEs, trait_encoded, 
                          method = "spearman", use = "p")  # Pearson 相关性

modTraitP <- corPvalueStudent(modTraitCor, nSamples = nrow(MEs_df))

my_title = "Module score and TNM grade relationships"

# 构建注释文本
text_mat <- paste(signif(modTraitCor, 2), "(", signif(modTraitP, 1), ")", sep = "")
dim(text_mat) = dim(modTraitCor)