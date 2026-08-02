setwd("C:/Users/26355/Desktop/xxx/niche breadth/")

library(spaa)
library(tidyverse)
library(readr)

# ===============================
# Step 1｜读入原始表格
# ===============================
asv_raw <- read.table("feature-table.tsv",
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE,
                      comment.char = "")

meta <- read_tsv("metadata.tsv") %>%
  rename(SampleID = `sample-id`) %>%
  mutate(
    treatment = factor(treatment, 
                       levels = c("Control", "PM", "A", "PMA", "Fe")),
    day = paste0("d", day)
  ) %>%
  mutate(day = factor(day, 
                      levels = c("d1", "d3", "d5", "d9", "d16", "d20")))

# ===============================
# Step 2｜筛选
# ===============================

otu <- asv_raw

# 👉 加在这里👇
otu_rel <- t(t(otu) / colSums(otu))

keep_asvs <- rowSums(otu_rel > 0.001) >= 5

otu <- otu[keep_asvs, ]

cat("ASV数量:", nrow(otu), "\n")

# ===============================
# Step 3｜计算每个OTU的niche width
# ===============================
# Spaa 要：行=样本，列=OTU
otu_t <- t(otu)

nb <- niche.width(otu_t, method = "levins")

# 注意：niche.width()返回的是data.frame，需要正确处理
nb_df <- data.frame(
  OTU = names(nb),
  niche_width = as.numeric(nb),  # 用as.numeric提取数值
  stringsAsFactors = FALSE
)

cat("\nniche_width计算完成，前5个OTU:\n")
print(head(nb_df))

# ===============================
# Step 4｜把OTU转到sample水平
# ===============================
# OTU表转长格式
otu_long <- as.data.frame(otu) %>%
  rownames_to_column("OTU") %>%
  pivot_longer(-OTU, names_to = "SampleID", values_to = "abundance")

# 合并niche width
otu_long <- left_join(otu_long, nb_df, by = "OTU")

# 合并metadata
otu_long <- left_join(otu_long, meta, by = "SampleID")

# 计算每个样本的mean_nb
# 注意：我们只计算在该样本中出现的OTU（abundance > 0）的niche_width平均值
nb_sample <- otu_long %>%
  filter(abundance > 0) %>%  # 只考虑出现的OTU
  group_by(SampleID, treatment) %>%
  summarise(
    mean_nb = mean(niche_width, na.rm = TRUE),
    n_otus = n(),  # 这个样本中出现的OTU数量
    .groups = "drop"
  )

cat("\n样本水平数据:\n")
cat("样本数量:", nrow(nb_sample), "\n")
cat("每个处理组的样本数量:\n")
print(table(nb_sample$treatment))

# 检查是否有样本没有OTU出现
if(any(is.na(nb_sample$mean_nb) | is.nan(nb_sample$mean_nb))) {
  cat("警告：有些样本的mean_nb是NA或NaN，将用0替换\n")
  nb_sample$mean_nb[is.na(nb_sample$mean_nb) | is.nan(nb_sample$mean_nb)] <- 0
}

# 检查每个样本的OTU数量
cat("\n每个样本中出现的OTU数量统计:\n")
print(summary(nb_sample$n_otus))

# ===============================
# Step 5｜ANOVA
# ===============================
anova_res <- aov(mean_nb ~ treatment, data = nb_sample)
anova_summary <- summary(anova_res)

# 提取统计值
stat_df <- as.data.frame(anova_summary[[1]])

df_between <- stat_df$Df[1]
df_within <- stat_df$Df[2]
f_value <- stat_df$`F value`[1]
p_value <- stat_df$`Pr(>F)`[1]

cat("\nANOVA结果:\n")
print(anova_summary)
cat("\n自由度: 组间 =", df_between, "，组内 =", df_within, 
    "，总自由度 =", df_between + df_within, "\n")
cat("总样本数 =", df_between + df_within + 1, "\n")

# 格式化p值显示
p_format <- ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value))
p_label <- paste0("ANOVA\nF(", df_between, ", ", df_within, ") = ",
                  round(f_value, 2),
                  "\np=", p_format)
# ===============================
# Step 6｜事后检验
# ===============================
if (p_value < 0.05) {
  library(emmeans)
  # 使用Tukey HSD多重比较
  tukey_res <- TukeyHSD(anova_res)
  cat("\nTukey HSD事后检验:\n")
  print(tukey_res)
  
  # 保存Tukey结果
  tukey_df <- as.data.frame(tukey_res$treatment)
  write.csv(tukey_df, "tukey_hsd_results_sample_level1.csv", row.names = TRUE)
  cat("\nTukey结果已保存到: tukey_hsd_results_sample_level1.csv\n")
}

# ===============================
# Step 7｜画图
# ===============================
# 设置颜色
group_colors <- c(
  "Control" = "#4FCB8B",  # 绿色
  "PM"      = "#4FC3E0",  # 蓝色
  "A"       = "#F4C96A",  # 黄色
  "PMA"     = "#E57373",  # 红色
  "Fe"      = "#B39DDB"   # 紫色
)

# 计算y轴范围
y_max <- max(nb_sample$mean_nb, na.rm = TRUE)
y_min <- min(nb_sample$mean_nb, na.rm = TRUE)

# 1. 首先计算每个处理组的中位数
median_values <- nb_sample %>%
  group_by(treatment) %>%
  summarise(median_nb = median(mean_nb, na.rm = TRUE))

# 2. 修改您的绘图代码 (p <- ggplot(...)部分)
p <- ggplot(nb_sample, aes(x = treatment, y = mean_nb)) +
  # 原始箱线图 (保持边框不变)
  geom_boxplot(aes(fill = treatment), 
               outlier.shape = NA, 
               coef = 0,
               alpha = 0.8,
               linewidth = 0.5, # 这是箱体边框的粗细
               fatten = NULL) + # 禁用默认的中位数线
  # ★ 使用geom_segment手动绘制中位线 (稳定可靠)
  geom_segment(data = median_values,
               aes(x = as.numeric(factor(treatment)) - 0.375, # 线段起点X (左移)
                   xend = as.numeric(factor(treatment)) + 0.375, # 线段终点X (右移)
                   y = median_nb,
                   yend = median_nb),
               color = "black", # 中位线颜色
               linewidth = 0.5, # ★ 中位线粗细，随意调整
               alpha = 1) + # 确保完全不透明
  # 3. 后续的散点、均值点等图层保持不变
  geom_jitter(color = "grey40", 
              width = 0.2, 
              size = 2, 
              alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", color = "black", size = 3) +
  # 将原代码行替换为：
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  
  # 填充颜色
  scale_fill_manual(values = group_colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 16, margin = margin(r = 10))
  ) +
  ylab("Niche breadth") +
  xlab("") +
  # 添加统计标注
  annotate("text", 
           x = Inf, y = Inf,
           label = p_label,
           hjust = 1.1, vjust = 1.5, size = 3.5) +
  # 可选：添加样本数量标签
  annotate("text", 
           x = 1:length(unique(nb_sample$treatment)), 
           y = y_min - 0.05 * (y_max - y_min),
           label = paste0("n=", as.numeric(table(nb_sample$treatment))),
           size = 3)

# 显示图形
print(p)

# 保存图形
ggsave("niche_breadth_plot_final_sample_level1.pdf", p, width = 7, height = 5)
ggsave("niche_breadth_plot_final_sample_level1.png", p, width = 7, height = 5, dpi = 300)

# 保存数据
write.csv(nb_sample, "sample_niche_breadth_final.csv", row.names = FALSE)

# ===============================
# Step 8｜统计摘要
# ===============================
group_stats <- nb_sample %>%
  group_by(treatment) %>%
  summarise(
    n_samples = n(),
    mean_niche_breadth = mean(mean_nb, na.rm = TRUE),
    sd_niche_breadth = sd(mean_nb, na.rm = TRUE),
    se_niche_breadth = sd_niche_breadth / sqrt(n_samples),
    median_niche_breadth = median(mean_nb, na.rm = TRUE),
    min_niche_breadth = min(mean_nb, na.rm = TRUE),
    max_niche_breadth = max(mean_nb, na.rm = TRUE)
  )

cat("\n各处理组生态位宽度统计（样本水平）:\n")
print(group_stats)
write.csv(group_stats, "group_statistics_final1.csv", row.names = FALSE)

# 完成
cat("\n分析完成！\n")
cat("图形已保存为: niche_breadth_plot_final_sample_level1.pdf/png\n")
cat("数据已保存为: sample_niche_breadth_final1.csv\n")