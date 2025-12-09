library(reshape2)
library(ggplot2)





# 创建数据框
data <- data.frame(
  Category = c("LC", "LT"),
  Mmul = c(0.021, 0.130),
  Sscrofa = c(0.979, 0.870)
)

# 将数据转换为长格式
data_long <- melt(data, id.vars = "Category")

# 绘制条形图
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/map/mRNA_mapping_barplot.pdf", height = 3, width = 4)
ggplot(data_long, aes(x = Category, y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 0.6, position = "dodge") +
  labs(title = "mRNA mapping rate",
       x = "GROUPS",
       y = "Values",
       fill = "Species") +
  scale_fill_manual(values = c("Mmul" = "#3B035F", "Sscrofa" = "#F8D33E")) +
  theme_classic()+
  theme(panel.grid.major = element_blank(),  # 去掉主要网格线
        panel.grid.minor = element_blank(),  # 去掉次要网格线
        axis.text.x = element_text(color = "black"),  # x 轴字体颜色
        axis.text.y = element_text(color = "black"),  # y 轴字体颜色
        axis.title.x = element_text(color = "black"),  # x 轴标题颜色
        axis.title.y = element_text(color = "black"),  # y 轴标题颜色
        axis.line = element_line(color = "black"),  # 添加坐标轴线
        axis.ticks.length = unit(0.1, "cm"))+
  scale_y_continuous(expand = c(0,0)) # 设置坐标轴刻度线长度
dev.off()




data <- data.frame(
  Category = c("LC", "LT"),
  Mmul = c(0.00061941, 0.12170364),
  Sscrofa = c(0.99938059, 0.87829636)
)

# 将数据转换为长格式
data_long <- melt(data, id.vars = "Category")

# 绘制条形图
pdf("/cluster/home/xyzhang_jh/projects/liver/analysis/wangrongrong/Xeno_pig_monkey/map/WGS_mapping_barplot.pdf", height = 3, width = 4)
ggplot(data_long, aes(x = Category, y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 0.6, position = "dodge") +
  labs(title = "WGS mapping rate",
       x = "GROUPS",
       y = "Values",
       fill = "Species") +
  scale_fill_manual(values = c("Mmul" = "#3B035F", "Sscrofa" = "#F8D33E")) +
  theme_classic()+
  theme(panel.grid.major = element_blank(),  # 去掉主要网格线
        panel.grid.minor = element_blank(),  # 去掉次要网格线
        axis.text.x = element_text(color = "black"),  # x 轴字体颜色
        axis.text.y = element_text(color = "black"),  # y 轴字体颜色
        axis.title.x = element_text(color = "black"),  # x 轴标题颜色
        axis.title.y = element_text(color = "black"),  # y 轴标题颜色
        axis.line = element_line(color = "black"),  # 添加坐标轴线
        axis.ticks.length = unit(0.1, "cm"))+
  scale_y_continuous(expand = c(0,0)) # 设置坐标轴刻度线长度
dev.off()