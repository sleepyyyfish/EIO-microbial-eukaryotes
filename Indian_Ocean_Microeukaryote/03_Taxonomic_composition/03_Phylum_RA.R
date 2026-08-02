setwd("C:/Users/26355/Desktop/xxx/relative abundance")

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
# Step 2｜清理 taxonomy（只做最必要的）
# ===============================
tax$Phylum <- gsub(".*__", "", tax$Phylum)
tax$Phylum <- trimws(tax$Phylum)

# NA 保留（不要乱改！）
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
  data.frame(OTU_ID = rownames(tax), Phylum = tax$Phylum),
  by = "OTU_ID"
)

otu_meta <- left_join(
  otu_tax,
  metadata,
  by = "sample-id"
)

# ===============================
# Step 5｜汇总到 Phylum（先按 replicate）
# ===============================
phylum_sum <- otu_meta %>%
  group_by(Phylum, day, treatment, replicate) %>%
  summarise(
    abundance = sum(abundance),
    .groups = "drop"
  )

# ===============================
# Step 6｜replicate 平均
# ===============================
phylum_mean <- phylum_sum %>%
  group_by(Phylum, day, treatment) %>%
  summarise(
    mean_abundance = mean(abundance),
    .groups = "drop"
  )



# ===============================
# Step 7｜选 Top11 已知门类（排除Unclassified）
# ===============================
top_phyla <- phylum_mean %>%
  group_by(Phylum) %>%
  summarise(total = sum(mean_abundance), .groups = "drop") %>%
  arrange(desc(total))

top11_known <- top_phyla %>%
  filter(Phylum != "Unclassified") %>%
  slice(1:11) %>%
  pull(Phylum)

# 计算占比
others_percent <- 1 - sum(top_phyla %>%
                            filter(Phylum %in% top11_known) %>%
                            pull(total)) / sum(top_phyla$total)

cat(sprintf("\nTop11已知门类占比: %.2f%%\n", (1-others_percent)*100))
cat(sprintf("Others占比（包含Unclassified）: %.2f%%\n", others_percent*100))

# ===============================
# Step 8｜Top11已知门类 + Others
# ===============================
phylum_plot <- phylum_mean %>%
  mutate(
    Phylum_group = ifelse(Phylum %in% top11_known, Phylum, "Others")
  ) %>%
  group_by(Phylum_group, day, treatment) %>%
  summarise(
    mean_abundance = sum(mean_abundance),
    .groups = "drop"
  )

# ===============================
# Step 9｜转相对丰度
# ===============================
phylum_plot_rel <- phylum_plot %>%
  group_by(day, treatment) %>%
  mutate(
    rel_abundance = mean_abundance / sum(mean_abundance)
  ) %>%
  ungroup()

# ===============================
# 修正：定义正确的堆积顺序
# ===============================
# 注意：ggplot2中堆积柱状图的顺序是从下到上，所以要反过来
# 我们希望自养在底部，混合在中间，异养在上面，Others在最顶部
# 所以定义顺序时，自养在最前面，混合其次，异养再次，Others在最后

phy_order_correct <- c(
  "Chlorophyta", "Diatomea",
  "Dinoflagellata","Cercozoa",  "Ochrophyta",
  "Arthropoda", "Bicosoecida", "MAST-2", "MAST-3", "Retaria",
  "Protalveolata",   # ✅ 寄生类
  "Others"
)

# 定义颜色（按照上面的顺序）
color_palette <- c(
  # Autotrophs
  "Chlorophyta" = "#447542",
  "Diatomea" = "#709020",
  # Mixotrophs
  "Dinoflagellata" = "#F4B842",
  "Cercozoa" = "#E49832",
  "Ochrophyta" = "#D47822",
  # Heterotrophs
  "Arthropoda" = "#F8CCD8",
  "Bicosoecida" = "#E8AAB8",
  "MAST-2" = "#D88A98",
  "MAST-3" = "#B86A78",
  "Retaria" = "#984A58",
  # Parasite（蓝色）
  "Protalveolata" = "#4A90b2",   # ✅ 蓝色
  # Others
  "Others" = "grey70"
)

# 定义因子顺序
phylum_plot_rel$Phylum_group <- factor(
  phylum_plot_rel$Phylum_group,
  levels = phy_order_correct
)

# ===============================
# Step 10｜处理因子顺序
# ===============================
phylum_plot_rel$day <- factor(
  phylum_plot_rel$day,
  levels = c("d1", "d3", "d5", "d9", "d16", "d20")
)

phylum_plot_rel$treatment <- factor(
  phylum_plot_rel$treatment,
  levels = c("Control", "PM", "A", "PMA", "Fe")
)

# ===============================
# Step 11｜画图
# ===============================
p <- ggplot(
  phylum_plot_rel,
  aes(
    x = treatment,
    y = rel_abundance,
    fill = Phylum_group
  )
) +
  # 修改1：去掉黑色边框和线条
  geom_bar(stat = "identity", width = 0.75) +
  facet_grid(
    ~ day,
    scales = "free_x",
    space = "free_x"
  ) +
  # 使用自定义颜色，并设置图例顺序
  scale_fill_manual(values = color_palette, 
                    breaks = phy_order_correct) +
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
    # 修改2：图例移到右侧
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.box.spacing = unit(0.3, "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    # 修改3：保留面板边框
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Phylum"
  ) +
  # 修改4：调整图例布局
  guides(fill = guide_legend(ncol = 1))

print(p)

# ===============================
# Step 12｜保存
# ===============================
ggsave("phylum_top11_ordered_relative.png", p, width = 12, height = 6, dpi = 300)
ggsave("phylum_top11_ordered_relative.pdf", p, width = 12, height = 6)

# ===============================
# 【修正】Step 13|计算并输出整体相对丰度汇总
# ===============================

# 1. 先计算出正确的「整体相对丰度」
# 思路：在整个数据集中，每个门的累计丰度占总丰度的比例
phylum_summary <- phylum_plot %>%                # 注意：这里用 phylum_plot（绝对丰度）
  group_by(Phylum_group) %>%
  summarise(total_abundance = sum(mean_abundance)) %>%
  mutate(overall_rel_abundance = total_abundance / sum(total_abundance))

# 2. 定义营养方式（和之前一致）
nutrition_groups <- list(
  Autotroph = c("Chlorophyta", "Diatomea"),
  Mixotroph = c("Cercozoa", "Dinoflagellata", "Ochrophyta"),
  Heterotroph = c("Arthropoda", "Bicosoecida", "MAST-2", "MAST-3", "Retaria"),
  Parasite = c("Protalveolata")   # 新增寄生类
)

# 3. 按你设定的顺序输出 Top11 + Others
cat("\n========== Top 11 Known Phyla + Others (Ordered) ==========\n")
cat("Phylum\t\tNutrition Mode\t\tOverall Relative Abundance\n")
cat("------------------------------------------------------------\n")

for(phy in phy_order_correct) {  # phy_order_correct 是你前面定义的颜色顺序
  summary_row <- phylum_summary %>%
    filter(Phylum_group == phy)
  
  if(nrow(summary_row) > 0) {
    rel_abund_val <- summary_row$overall_rel_abundance[1]
    rel_abund_str <- paste0(format(round(rel_abund_val * 100, 2), nsmall = 2), "%")
    
    # 判断营养方式
    nutrition_mode <- ifelse(phy %in% nutrition_groups$Autotroph, "Autotroph",
                             ifelse(phy %in% nutrition_groups$Mixotroph, "Mixotroph",
                                    ifelse(phy %in% nutrition_groups$Heterotroph, "Heterotroph",
                                           ifelse(phy %in% nutrition_groups$Parasite, "Parasite", "Other/Unknown"))))
    
    cat(sprintf("%-12s\t%-15s\t%20s\n", 
                ifelse(nchar(phy) > 10, paste0(substr(phy, 1, 8), ".."), phy),
                nutrition_mode,
                rel_abund_str))
  }
}
# ==================================================
# 【修正】Step 14｜导出：每个处理×时间的相对丰度宽表
# ==================================================

# 1️⃣ 锁定 day / treatment 顺序（与你画图完全一致）
phylum_plot_rel$day <- factor(
  phylum_plot_rel$day,
  levels = c("d1", "d3", "d5", "d9", "d16", "d20")
)

phylum_plot_rel$treatment <- factor(
  phylum_plot_rel$treatment,
  levels = c("Control", "PM", "A", "PMA", "Fe")
)

# 2️⃣ 按 day × treatment 排序（关键！）
phylum_plot_rel <- phylum_plot_rel %>%
  arrange(day, treatment)

# 3️⃣ 构建分组列
phylum_export <- phylum_plot_rel %>%
  mutate(group_col = paste(day, treatment, sep = "_")) %>%
  select(Phylum_group, group_col, rel_abundance)

# 4️⃣ 转宽表（列顺序现在已固定）
wide_for_prism <- phylum_export %>%
  pivot_wider(
    names_from = group_col,
    values_from = rel_abundance,
    values_fill = 0
  )

# 5️⃣ 行顺序 = 图例顺序（phy_order_correct）
wide_for_prism <- wide_for_prism %>%
  mutate(Phylum_group = factor(Phylum_group, levels = phy_order_correct)) %>%
  arrange(Phylum_group)

# 6️⃣ 保存
write.csv(
  wide_for_prism,
  file = "phylum_top11_relative_abundance_for_Prism.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("\n✅ 已导出：phylum_top11_relative_abundance_for_Prism.csv\n")
cat("  行顺序 = 图例顺序；列顺序 = d1→d20 × Control→Fe\n")