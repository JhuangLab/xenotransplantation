library(reshape2)
library(ggplot2)


gene_names <- c("FN1","FGG","ITGA5","PECAM1","UBA3","GRN","CISD1","UBA2","UBE4A","FGA","MPP1","ITGB1","GSN",
                "MAP2K1","RAP1A","RALA","ARF6","MAPRE1","RTN4","MACF1","MAP2K3","PKN1","UBE2K","HSPA8","PSMA2","PSMA5")

count_mmul <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm_mmul.txt")
# 去掉有0的行
count_mmul <- count_mmul[rownames(count_mmul) %in% gene_names,]
# 去掉全是 0 的行
#count_total <- count_total[rowSums(count_total != 0) > 0, ]

count_mmul <- count_mmul * 1e6
count_mmul <- log10(count_mmul+1)
count_mmul <- as.data.frame(count_mmul)

# 计算每列的平均值
averages <- colMeans(count_mmul)  # 只计算LC和LT列的平均值

# 创建新的数据框来存储结果
result_mmul <- data.frame(
  LC = averages[1],
  LT = averages[2]
)
rownames(result_mmul) <- c("mmul")
result_mmul <- t(result_mmul) %>% as.data.frame()
result_mmul$group <- rownames(result_mmul)
# 将数据从宽格式转换为长格式
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/mmul_up_gene_barplot.pdf",height = 3, width = 3)
ggplot(result_mmul, aes(x = group, y = mmul, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, position = "dodge") +
  labs(title = "Mmul_gene_expression",
       x = "Group",
       y = "Values",
       fill = "Category")+
  scale_fill_manual(values = c("LC" = "#3B035F", "LT" = "#F8D33E")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),  # 去掉主要网格线
        panel.grid.minor = element_blank(),  # 去掉次要网格线
        axis.text.x = element_text(color = "black"),  # x 轴字体颜色
        axis.text.y = element_text(color = "black"),  # y 轴字体颜色
        axis.title.x = element_text(color = "black"),  # x 轴标题颜色
        axis.title.y = element_text(color = "black"),  # y 轴标题颜色
        axis.line = element_line(color = "black"),  # 添加坐标轴线
        axis.ticks.length = unit(0.1, "cm")) +
  scale_y_continuous(expand = c(0, 0))+# 设置坐标轴刻度线长度
  NoLegend()
dev.off()

##sscrofa-----------------------------------------------------------------------------
count_sscrofa <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm_sscrofa.txt")

# 去掉有0的行
count_sscrofa <- count_sscrofa[rownames(count_sscrofa) %in% gene_names,]
# 去掉全是 0 的行
#count_total <- count_total[rowSums(count_total != 0) > 0, ]

count_sscrofa <- count_sscrofa * 1e6
count_sscrofa <- log10(count_sscrofa+1)
count_sscrofa <- as.data.frame(count_sscrofa)

# 计算每列的平均值
averages <- colMeans(count_sscrofa)  # 只计算LC和LT列的平均值

# 创建新的数据框来存储结果
result_sscrofa <- data.frame(
  LC = averages[1],
  LT = averages[2]
)
rownames(result_sscrofa) <- c("sscrofa")
result_sscrofa <- t(result_sscrofa) %>% as.data.frame()
result_sscrofa$group <- rownames(result_sscrofa)
# 将数据从宽格式转换为长格式
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/sscrofa_up_gene_barplot.pdf",height = 3, width = 3)
ggplot(result_sscrofa, aes(x = group, y = sscrofa, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, position = "dodge") +
  labs(title = "Sscrofa_gene_expression",
       x = "Group",
       y = "Values",
       fill = "Category")+
  scale_fill_manual(values = c("LC" = "#3B035F", "LT" = "#F8D33E")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),  # 去掉主要网格线
        panel.grid.minor = element_blank(),  # 去掉次要网格线
        axis.text.x = element_text(color = "black"),  # x 轴字体颜色
        axis.text.y = element_text(color = "black"),  # y 轴字体颜色
        axis.title.x = element_text(color = "black"),  # x 轴标题颜色
        axis.title.y = element_text(color = "black"),  # y 轴标题颜色
        axis.line = element_line(color = "black"),  # 添加坐标轴线
        axis.ticks.length = unit(0.1, "cm")) +
  scale_y_continuous(expand = c(0, 0))+# 设置坐标轴刻度线长度
  NoLegend()
dev.off()




##v2_2
all_exp <- read.delim("/cluster/home/jhuang/projects/liver/data/zhangwei/pig/protein/report/report-liver/data/3.comparison/PIT_CON-IT/PIT_CON-IT.all.xls")

all_exp <- all_exp[,c("gene_name","CON.IT","PIT")]
all_exp <- all_exp[!grepl("-", all_exp$gene_name), ] %>% na.omit()
all_exp <- all_exp[!duplicated(all_exp$gene_name),]
rownames(all_exp) <- NULL
all_exp <- all_exp %>% column_to_rownames(.,"gene_name")

column_sums <- colSums(all_exp)
normalized_protein <- t(t(all_exp) / column_sums) * 1e6  # 将标准化后的值乘以 1,000,000
normalized_protein <- log10(normalized_protein+1)
normalized_protein <- as.data.frame(normalized_protein)
normalized_protein <- normalized_protein[rownames(normalized_protein) %in% gene_names,]

averages <- colMeans(normalized_protein)  # 只计算LC和LT列的平均值

# 创建新的数据框来存储结果
result_protein <- data.frame(
  CON.IT = averages[1],
  PIT = averages[2]
)
rownames(result_protein) <- c("protein")
result_protein <- t(result_protein) %>% as.data.frame()
result_protein$group <- rownames(result_protein)
# 将数据从宽格式转换为长格式
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/protein_up_gene_barplot.pdf",height = 3, width = 3)
ggplot(result_protein, aes(x = group, y = protein, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, position = "dodge") +
  labs(title = "Protein_level",
       x = "Group",
       y = "Values",
       fill = "Category")+
  scale_fill_manual(values = c("CON.IT" = "#3B035F", "PIT" = "#F8D33E")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),  # 去掉主要网格线
        panel.grid.minor = element_blank(),  # 去掉次要网格线
        axis.text.x = element_text(color = "black"),  # x 轴字体颜色
        axis.text.y = element_text(color = "black"),  # y 轴字体颜色
        axis.title.x = element_text(color = "black"),  # x 轴标题颜色
        axis.title.y = element_text(color = "black"),  # y 轴标题颜色
        axis.line = element_line(color = "black"),  # 添加坐标轴线
        axis.ticks.length = unit(0.1, "cm")) +
  scale_y_continuous(expand = c(0, 0))+# 设置坐标轴刻度线长度
  NoLegend()
dev.off()

#down------------------------------------------------------------------------------------------------------
gene_names <- c("SULT1E1","SULT2A1","FTL","FMO1","ADH4","AOX1","CTSB","CTSS","PTCD3","MRPL50","HBB","SLC4A1","SLC2A2",
                "GYS2","GULO","TDH","SORD","SEPSECS","ETNPPL","PHGDH","TAT","TDO2","PCK1","FBP1","SDS","SDSL","ACACA","PPAT","DPYD")

count_mmul <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm_mmul.txt")
# 去掉有0的行
count_mmul <- count_mmul[rownames(count_mmul) %in% gene_names,]
# 去掉全是 0 的行
#count_total <- count_total[rowSums(count_total != 0) > 0, ]

count_mmul <- count_mmul * 1e6
count_mmul <- log10(count_mmul+1)
count_mmul <- as.data.frame(count_mmul)

# 计算每列的平均值
averages <- colMeans(count_mmul)  # 只计算LC和LT列的平均值

# 创建新的数据框来存储结果
result_mmul <- data.frame(
  LC = averages[1],
  LT = averages[2]
)
rownames(result_mmul) <- c("mmul")
result_mmul <- t(result_mmul) %>% as.data.frame()
result_mmul$group <- rownames(result_mmul)
# 将数据从宽格式转换为长格式
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/mmul_down_gene_barplot.pdf",height = 3, width = 3)
ggplot(result_mmul, aes(x = group, y = mmul, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, position = "dodge") +
  labs(title = "Mmul_gene_expression",
       x = "Group",
       y = "Values",
       fill = "Category")+
  scale_fill_manual(values = c("LC" = "#3B035F", "LT" = "#F8D33E")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),  # 去掉主要网格线
        panel.grid.minor = element_blank(),  # 去掉次要网格线
        axis.text.x = element_text(color = "black"),  # x 轴字体颜色
        axis.text.y = element_text(color = "black"),  # y 轴字体颜色
        axis.title.x = element_text(color = "black"),  # x 轴标题颜色
        axis.title.y = element_text(color = "black"),  # y 轴标题颜色
        axis.line = element_line(color = "black"),  # 添加坐标轴线
        axis.ticks.length = unit(0.1, "cm")) +
  scale_y_continuous(expand = c(0, 0))+# 设置坐标轴刻度线长度
  NoLegend()
dev.off()

##sscrofa-----------------------------------------------------------------------------
count_sscrofa <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts_fpkm_sscrofa.txt")

# 去掉有0的行
count_sscrofa <- count_sscrofa[rownames(count_sscrofa) %in% gene_names,]
# 去掉全是 0 的行
#count_total <- count_total[rowSums(count_total != 0) > 0, ]

count_sscrofa <- count_sscrofa * 1e6
count_sscrofa <- log10(count_sscrofa+1)
count_sscrofa <- as.data.frame(count_sscrofa)

# 计算每列的平均值
averages <- colMeans(count_sscrofa)  # 只计算LC和LT列的平均值

# 创建新的数据框来存储结果
result_sscrofa <- data.frame(
  LC = averages[1],
  LT = averages[2]
)
rownames(result_sscrofa) <- c("sscrofa")
result_sscrofa <- t(result_sscrofa) %>% as.data.frame()
result_sscrofa$group <- rownames(result_sscrofa)
# 将数据从宽格式转换为长格式
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/sscrofa_down_gene_barplot.pdf",height = 3, width = 3)
ggplot(result_sscrofa, aes(x = group, y = sscrofa, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, position = "dodge") +
  labs(title = "Sscrofa_gene_expression",
       x = "Group",
       y = "Values",
       fill = "Category")+
  scale_fill_manual(values = c("LC" = "#3B035F", "LT" = "#F8D33E")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),  # 去掉主要网格线
        panel.grid.minor = element_blank(),  # 去掉次要网格线
        axis.text.x = element_text(color = "black"),  # x 轴字体颜色
        axis.text.y = element_text(color = "black"),  # y 轴字体颜色
        axis.title.x = element_text(color = "black"),  # x 轴标题颜色
        axis.title.y = element_text(color = "black"),  # y 轴标题颜色
        axis.line = element_line(color = "black"),  # 添加坐标轴线
        axis.ticks.length = unit(0.1, "cm")) +
  scale_y_continuous(expand = c(0, 0))+# 设置坐标轴刻度线长度
  NoLegend()
dev.off()




##v2_2
all_exp <- read.delim("/cluster/home/jhuang/projects/liver/data/zhangwei/pig/protein/report/report-liver/data/3.comparison/PIT_CON-IT/PIT_CON-IT.all.xls")

all_exp <- all_exp[,c("gene_name","CON.IT","PIT")]
all_exp <- all_exp[!grepl("-", all_exp$gene_name), ] %>% na.omit()
all_exp <- all_exp[!duplicated(all_exp$gene_name),]
rownames(all_exp) <- NULL
all_exp <- all_exp %>% column_to_rownames(.,"gene_name")

column_sums <- colSums(all_exp)
normalized_protein <- t(t(all_exp) / column_sums) * 1e6  # 将标准化后的值乘以 1,000,000
normalized_protein <- log10(normalized_protein+1)
normalized_protein <- as.data.frame(normalized_protein)
normalized_protein <- normalized_protein[rownames(normalized_protein) %in% gene_names,]

averages <- colMeans(normalized_protein)  # 只计算LC和LT列的平均值

# 创建新的数据框来存储结果
result_protein <- data.frame(
  CON.IT = averages[1],
  PIT = averages[2]
)
rownames(result_protein) <- c("protein")
result_protein <- t(result_protein) %>% as.data.frame()
result_protein$group <- rownames(result_protein)
# 将数据从宽格式转换为长格式
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/venn/protein_down_gene_barplot.pdf",height = 3, width = 3)
ggplot(result_protein, aes(x = group, y = protein, fill = group)) +
  geom_bar(stat = "identity", width = 0.6, position = "dodge") +
  labs(title = "Protein_level",
       x = "Group",
       y = "Values",
       fill = "Category")+
  scale_fill_manual(values = c("CON.IT" = "#3B035F", "PIT" = "#F8D33E")) +
  theme_classic() +
  theme(panel.grid.major = element_blank(),  # 去掉主要网格线
        panel.grid.minor = element_blank(),  # 去掉次要网格线
        axis.text.x = element_text(color = "black"),  # x 轴字体颜色
        axis.text.y = element_text(color = "black"),  # y 轴字体颜色
        axis.title.x = element_text(color = "black"),  # x 轴标题颜色
        axis.title.y = element_text(color = "black"),  # y 轴标题颜色
        axis.line = element_line(color = "black"),  # 添加坐标轴线
        axis.ticks.length = unit(0.1, "cm")) +
  scale_y_continuous(expand = c(0, 0))+# 设置坐标轴刻度线长度
  NoLegend()
dev.off()
