.libPaths(new = c(.libPaths(), "/cluster/home/danyang_jh/sbin/R/R-4.3.0", "/cluster/home/danyang_jh/sbin/R/R-4.2.1","/cluster/home/xyzhang_jh/sbin/R"))
pkgs <- c("fs", "futile.logger", "configr", "stringr", "ggpubr", "ggthemes", "QFeatures", 
          "jhtools", "glue", "ggsci", "patchwork", "tidyverse", "dplyr", "MSnbase", "msqrob2", 
          "enrichplot", "org.Hs.eg.db", "parallel", "jhuanglabRNAseq", "jhuanglabGO", "limma", 
          "readxl", "writexl", "ComplexHeatmap", "circlize", "clusterProfiler", "MsCoreUtils", 
          "missForest", "paletteer", "factoextra", "FactoMineR", "openxlsx","gdata","ggrepel","edgeR")  
for (pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}


raw_counts <- read.table(file = "/cluster/home/yliang_jh/projects/mRNA/xenograft_zhangwei/results/salmon_quant/counts.txt")

#拆分物种
raw_counts$Species <- sapply(str_split(rownames(raw_counts), "_"), `[`, 1)

#计算所有物种的总和
species_LC_sum <- sum(raw_counts$LC)
species_LT_sum <- sum(raw_counts$LT)

#计算单一物种的总和
Mmul_table <- raw_counts[raw_counts$Species == "Mmul", ]
Mmul_LC_ratio <- sum(Mmul_table$LC)/species_LC_sum
Mmul_LT_ratio <- sum(Mmul_table$LT)/species_LT_sum


Sscrofa_table <- raw_counts[raw_counts$Species == "Sscrofa", ]
Sscrofa_LC_ratio <- sum(Sscrofa_table$LC)/species_LC_sum
Sscrofa_LT_ratio <- sum(Sscrofa_table$LT)/species_LT_sum

result <- data.frame(
  LC = c(Mmul_LC_ratio,Sscrofa_LC_ratio),
  LT = c(Mmul_LT_ratio,Sscrofa_LT_ratio)
)

rownames(result) <- c("Mmul", "Sscrofa")
result$species <- rownames(result)
write.csv(result, file = "/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/mRNA/mapping_ratio.csv")

table1_output <- table1(~ LC+LT | species, data = result, overall = FALSE)
print(table1_output)