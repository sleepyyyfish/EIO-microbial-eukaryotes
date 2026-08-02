setwd("C:/Users/26355/Desktop/毕设111/aaa重来一遍/relative abundance")

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ===============================
# Step 1｜读入数据
# ===============================
otu <- read.table("feature_table.tsv",
                  sep = "\t",
                  header = TRUE,
                  row.names = 1,
                  check.names = FALSE)

tax <- read.table("taxonomy.tsv",
                  sep = "\t",
                  header = TRUE,
                  row.names = 1,
                  fill = TRUE)

metadata <- read.table("metadata_final.tsv",
                       sep = "\t",
                       header = TRUE,
                       stringsAsFactors = FALSE,
                       check.names = FALSE)

# ===============================
# Step 2｜清理 taxonomy（Class级别）
# ===============================
tax$Class <- gsub(".*__", "", tax$Class)
tax$Class <- trimws(tax$Class)

# NA 和空值处理
tax$Class[is.na(tax$Class) | tax$Class == ""] <- "Unclassified"

# 也清理Phylum列，用于营养方式分类
tax$Phylum <- gsub(".*__", "", tax$Phylum)
tax$Phylum[is.na(tax$Phylum) | tax$Phylum == ""] <- "Unclassified"

# ===============================
# Step 3｜OTU → 长表
# ===============================
otu$OTU_ID <- rownames(otu)

otu_long <- pivot_longer(
  otu,
  cols = -OTU_ID,
  names_to = "sample-id",
  values_to = "abundance"
)

# ===============================
# Step 4｜合并 taxonomy + metadata
# ===============================
otu_tax <- left_join(
  otu_long,
  data.frame(OTU_ID = rownames(tax), 
             Class = tax$Class,
             Phylum = tax$Phylum),  # 保留Phylum用于营养方式分类
  by = "OTU_ID"
)

otu_meta <- left_join(
  otu_tax,
  metadata,
  by = "sample-id"
)

# ===============================
# Step 5｜汇总到 Class（先按 replicate）
# ===============================
class_sum <- otu_meta %>%
  group_by(Class, day, treatment, replicate) %>%
  summarise(
    abundance = sum(abundance),
    .groups = "drop"
  )

# ===============================
# Step 6｜replicate 平均
# ===============================
class_mean <- class_sum %>%
  group_by(Class, day, treatment) %>%
  summarise(
    mean_abundance = mean(abundance),
    .groups = "drop"
  )

# ===============================
# Step 7｜选 Top15 已知类别（排除Unclassified）
# ===============================
top_classes <- class_mean %>%
  group_by(Class) %>%
  summarise(total = sum(mean_abundance), .groups = "drop") %>%
  arrange(desc(total))

top15_known <- top_classes %>%
  filter(Class != "Unclassified") %>%
  slice(1:15) %>%
  pull(Class)

# 计算占比
others_percent <- 1 - sum(top_classes %>%
                            filter(Class %in% top15_known) %>%
                            pull(total)) / sum(top_classes$total)

cat(sprintf("\nTop15已知Class占比: %.2f%%\n", (1-others_percent)*100))
cat(sprintf("Others占比（包含Unclassified）: %.2f%%\n", others_percent*100))

# ===============================
# Step 8｜Top15已知类别 + Others
# ===============================
class_plot <- class_mean %>%
  mutate(
    Class_group = ifelse(Class %in% top15_known, Class, "Others")
  ) %>%
  group_by(Class_group, day, treatment) %>%
  summarise(
    mean_abundance = sum(mean_abundance),
    .groups = "drop"
  )

# ===============================
# Step 9｜转相对丰度
# ===============================
class_plot_rel <- class_plot %>%
  group_by(day, treatment) %>%
  mutate(
    rel_abundance = mean_abundance / sum(mean_abundance)
  ) %>%
  ungroup()

# ===============================
# Step 10｜定义营养方式和Class顺序（按你的要求）
# ===============================

# 创建Class到Phylum的映射
class_to_phylum <- otu_meta %>%
  distinct(Class, Phylum) %>%
  filter(!is.na(Class), !is.na(Phylum))

# 定义营养方式分类（基于Phylum）
# 特别注意：将Syndiniales（Protalveolata）归为Parasite
nutrition_groups <- list(
  Autotroph = c("Chlorophyta", "Diatomea","Prymnesiophyceae"),
  Mixotroph = c("Cercozoa", "Dinoflagellata", "Ochrophyta"),  # 移除了Protalveolata
  Heterotroph = c("Arthropoda", "Bicosoecida","Centrohelida","MAST-2", "MAST-3", "Retaria","SA1-3C06"),
  Parasite = c("Protalveolata")  # 新增寄生类
)

# 为每个Class分配营养方式
class_nutrition <- class_to_phylum %>%
  mutate(
    nutrition = case_when(
      Phylum %in% nutrition_groups$Autotroph ~ "Autotroph",
      Phylum %in% nutrition_groups$Mixotroph ~ "Mixotroph",
      Phylum %in% nutrition_groups$Heterotroph ~ "Heterotroph",
      Phylum %in% nutrition_groups$Parasite ~ "Parasite",  # 新增寄生类
      TRUE ~ "Other/Unclassified"
    )
  )

# ===============================
# 关键修改：按你指定的顺序排列Class
# ===============================
# 你指定的顺序：
# 自养：Clade_VII ，Bacillariophyceae, Mediophyceae，Prymnesiophyceae
# 混合：Dinophyceae, Chlorarachniophyta，Dictyochophyceae
# 异养：Maxillopoda, Bicosoecida, MAST-2, MAST-3C, RAD_B, Centrohelida, SA1-3C06
# 寄生：Syndiniales
# Others

desired_order <- c(
  # 自养（绿色系）
  "Clade_VII", "Bacillariophyceae", "Mediophyceae", "Prymnesiophyceae",
  # 混合（黄色/橙色系）
  "Dinophyceae", "Chlorarachniophyta", "Dictyochophyceae",
  # 异养（红色系）
  "Maxillopoda", "Bicosoecida", "MAST-2", "MAST-3C", "RAD_B", "Centrohelida", "SA1-3C06",
  # 寄生（蓝色系）
  "Syndiniales",
  # Others
  "Others"
)

# 只保留Top15中存在的Class，并保持指定顺序
class_order_correct <- desired_order[desired_order %in% top15_known]

# 如果有Top15中的Class不在指定顺序中，添加到Others之前
missing_classes <- setdiff(top15_known, desired_order)
if(length(missing_classes) > 0) {
  class_order_correct <- c(class_order_correct, missing_classes, "Others")
} else {
  class_order_correct <- c(class_order_correct, "Others")
}

# ===============================
# Step 11｜定义颜色方案（按营养方式）
# ===============================
# 为不同的营养方式定义颜色系
autotroph_colors <- c("#1B9E77", "#81C784", "#A5D6A7", "#C5E1A5")  # 绿色系
mixotroph_colors <- c("#FFEB3B", "#FBC02D", "#F9A825", "#F68D00")  # 黄色/橙色系
heterotroph_colors <- c("#FFD5C5", "#FF9C85", "#DAABCB", "#C868A1","#D33F6A","#E29EAB","#D66982")  # 红色系
parasite_colors <- c("#7AC7E2")  # 蓝色系（寄生）

# 创建颜色向量
color_palette_class <- c()

# 为每个Class分配颜色
for(class in class_order_correct) {
  if(class == "Others") {
    color_palette_class[class] <- "grey90"
    next
  }
  
  nutrition <- class_nutrition$nutrition[class_nutrition$Class == class]
  
  if(nutrition == "Autotroph" && length(autotroph_colors) > 0) {
    color_palette_class[class] <- autotroph_colors[1]
    autotroph_colors <- autotroph_colors[-1]
  } else if(nutrition == "Mixotroph" && length(mixotroph_colors) > 0) {
    color_palette_class[class] <- mixotroph_colors[1]
    mixotroph_colors <- mixotroph_colors[-1]
  } else if(nutrition == "Heterotroph" && length(heterotroph_colors) > 0) {
    color_palette_class[class] <- heterotroph_colors[1]
    heterotroph_colors <- heterotroph_colors[-1]
  } else if(nutrition == "Parasite" && length(parasite_colors) > 0) {
    color_palette_class[class] <- parasite_colors[1]  # 蓝色
    parasite_colors <- parasite_colors[-1]
  } else {
    # 如果颜色用完了，使用默认颜色
    color_palette_class[class] <- "#CCCCCC"
  }
}

# ===============================
# Step 12｜处理因子顺序
# ===============================
class_plot_rel$Class_group <- factor(
  class_plot_rel$Class_group,
  levels = class_order_correct
)

class_plot_rel$day <- factor(
  class_plot_rel$day,
  levels = c("d1", "d3", "d5", "d9", "d16", "d20")
)

class_plot_rel$treatment <- factor(
  class_plot_rel$treatment,
  levels = c("Control", "PM", "A", "PMA", "Fe")
)

# ===============================
# Step 13｜画图
# ===============================
p <- ggplot(
  class_plot_rel,
  aes(
    x = treatment,
    y = rel_abundance,
    fill = Class_group
  )
) +
  geom_bar(stat = "identity", width = 0.75) +
  facet_grid(
    ~ day,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_fill_manual(values = color_palette_class) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_bw() +
  theme(
    panel.spacing.x = unit(0.8, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text.y = element_text(size = 15),
    axis.title = element_text(face = "bold", size = 15),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 12),
    legend.box.spacing = unit(0.3, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Class"
  ) +
  guides(fill = guide_legend(ncol = 1))

print(p)

# ===============================
# Step 14｜保存
# ===============================
ggsave("class_top15_ordered_relative.png", p, width = 12, height = 6, dpi = 300)
ggsave("class_top15_ordered_relative.pdf", p, width = 12, height = 6)

# ===============================
# Step 15｜输出Top15已知Class汇总
# ===============================
cat("\n========== Top 15 Known Classes (Custom Order) ==========\n")
cat("Class\t\tPhylum\t\t\tNutrition Mode\t\tOverall Relative Abundance\n")
cat("-------------------------------------------------------------------------\n")

# 计算每个Class的整体相对丰度
class_summary <- class_mean %>%
  group_by(Class) %>%
  summarise(total_abundance = sum(mean_abundance), .groups = "drop") %>%
  mutate(
    overall_rel_abundance = total_abundance / sum(total_abundance)
  ) %>%
  left_join(class_nutrition, by = "Class")

# 按自定义顺序输出
for(class_name in class_order_correct) {
  if(class_name != "Others") {
    summary_row <- class_summary %>%
      filter(Class == class_name)
    
    if(nrow(summary_row) > 0) {
      rel_abund <- paste0(format(round(summary_row$overall_rel_abundance[1] * 100, 2), nsmall = 2), "%")
      phylum_name <- summary_row$Phylum[1]
      nutrition_mode <- summary_row$nutrition[1]
      
      cat(sprintf("%-15s\t%-15s\t%-15s\t%20s\n", 
                  ifelse(nchar(class_name) > 15, paste0(substr(class_name, 1, 13), ".."), class_name),
                  ifelse(nchar(phylum_name) > 15, paste0(substr(phylum_name, 1, 13), ".."), phylum_name),
                  nutrition_mode,
                  rel_abund))
    }
  }
}

# 添加Others信息
others_abund <- class_summary %>%
  filter(!Class %in% top15_known) %>%
  summarise(total = sum(overall_rel_abundance))

cat("\nOthers (includes Unclassified and other minor groups):\n")
cat(sprintf("Overall Relative Abundance: %s\n", 
            paste0(format(round(others_abund$total * 100, 2), nsmall = 2), "%")))

# 输出营养方式统计
cat("\n========== Nutrition Mode Summary ==========\n")
nutrition_summary <- class_summary %>%
  filter(Class %in% top15_known) %>%
  group_by(nutrition) %>%
  summarise(
    count = n(),
    total_rel_abundance = sum(overall_rel_abundance) * 100
  )

for(i in 1:nrow(nutrition_summary)) {
  cat(sprintf("%s: %d classes, %.2f%%\n", 
              nutrition_summary$nutrition[i],
              nutrition_summary$count[i],
              nutrition_summary$total_rel_abundance[i]))
}
# ==================================================
# 【修正】Step 16｜导出：Class 级别相对丰度宽表
# ==================================================

# 1️⃣ 锁定 day / treatment 顺序
class_plot_rel$day <- factor(
  class_plot_rel$day,
  levels = c("d1", "d3", "d5", "d9", "d16", "d20")
)

class_plot_rel$treatment <- factor(
  class_plot_rel$treatment,
  levels = c("Control", "PM", "A", "PMA", "Fe")
)

# 2️⃣ 按 day × treatment 排序（关键）
class_plot_rel <- class_plot_rel %>%
  arrange(day, treatment)

# 3️⃣ 构建分组列
class_export <- class_plot_rel %>%
  mutate(group_col = paste(day, treatment, sep = "_")) %>%
  select(Class_group, group_col, rel_abundance)

# 4️⃣ 转宽表
wide_class_for_prism <- class_export %>%
  pivot_wider(
    names_from = group_col,
    values_from = rel_abundance,
    values_fill = 0
  )

# 5️⃣ 行顺序 = 图例顺序（class_order_correct）
wide_class_for_prism <- wide_class_for_prism %>%
  mutate(Class_group = factor(Class_group, levels = class_order_correct)) %>%
  arrange(Class_group)

# 6️⃣ 保存
write.csv(
  wide_class_for_prism,
  file = "class_top15_relative_abundance_for_Prism.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("\n✅ 已导出：class_top15_relative_abundance_for_Prism.csv\n")
cat("  行顺序 = 图例顺序；列顺序 = d1→d20 × Control→Fe\n")