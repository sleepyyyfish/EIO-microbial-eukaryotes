# ======================
# 0. 设置路径 + 加载包
# ======================
setwd("C:/Users/26355/Desktop/xxx/SIMPER")

library(vegan)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(stringr)
library(scales)

# ======================
# 1. 读入数据
# ======================

meta <- read_tsv("metadata_final.tsv")
tax <- read_tsv("taxonomy.tsv")

# ======================
# 2. 整理feature table
# ======================

otu <- read_tsv("feature_table_RA.tsv")
otu_ids <- otu$OTU_ID
rownames(otu) <- otu$OTU_ID
otu$OTU_ID <- NULL
otu_t <- t(otu)
otu_df <- as.data.frame(otu_t)
colnames(otu_df) <- otu_ids
otu <- otu_df
otu <- otu[meta$`sample-id`, ]

# ======================
# 3. 定义所有比较组
# ======================

comparisons <- list(
  Control_PM = c("Control", "PM"),
  Control_A = c("Control", "A"),
  Control_PMA = c("Control", "PMA"),
  Control_Fe = c("Control", "Fe"),
  PM_PMA = c("PM", "PMA"),
  A_PMA = c("A", "PMA")
)

group_levels <- c("Control", "PM", "A", "PMA", "Fe")
group_colors <- c(
  "Control" = "#82C262",
  "PM" = "#5BA2CC",
  "A" = "#F4D166",
  "PMA" = "#F17060",
  "Fe" = "#AC6A9F"
)

point_color <- "#BD1100"

# ======================
# 4. 修正的SIMPER分析函数
# ======================

run_corrected_simper <- function(comp_name, group1, group2, threshold = 0.02) {
  cat("\n=== 分析:", comp_name, "(", group1, "vs", group2, ") ===\n")

  # 筛选数据
  groups <- c(group1, group2)
  meta_sub <- meta %>% filter(treatment %in% groups)
  otu_sub <- otu[rownames(otu) %in% meta_sub$`sample-id`, ]
  meta_sub <- meta_sub[match(rownames(otu_sub), meta_sub$`sample-id`), ]

  # 运行SIMPER
  sim <- simper(otu_sub, as.character(meta_sub$treatment))

  # 提取比较结果
  simper_name <- paste0(group1, "_", group2)
  res <- sim[[simper_name]]
  df <- as.data.frame(res)
  df$OTU_ID <- rownames(df)

  # 重新计算关键指标，避免混淆
  # SIMPER原始输出中已经有average、contri、cumsum列
  # 但我们重新计算以确保准确性

  # 计算总贡献值
  total_average <- sum(df$average, na.rm = TRUE)

  # 重新计算贡献比例
  df$contribution_prop <- df$average / total_average

  # 按贡献值排序
  df <- df %>% arrange(desc(average))

  # 重新计算累计贡献比例
  df$cumulative_prop <- cumsum(df$contribution_prop)

  # 筛选贡献>阈值的OTU
  df_top <- df %>% filter(contribution_prop > threshold)

  if (nrow(df_top) == 0) {
    cat("  警告: 没有OTU满足筛选条件，使用前10个\n")
    df_top <- df %>% head(10)
  }

  # 合并分类学信息
  df_top <- left_join(df_top, tax, by = "OTU_ID")

  # 生成标签
  df_top <- df_top %>%
    mutate(
      # 清理分类学信息
      across(c(Phylum, Class, Order, Family, Genus, Species), 
             ~ifelse(is.na(.) | . == "", NA, .)),

      # 去除前缀
      Phylum_clean = str_replace(Phylum, "^[^_]*__", ""),
      Class_clean = str_replace(Class, "^[^_]*__", ""),
      Order_clean = str_replace(Order, "^[^_]*__", ""),
      Family_clean = str_replace(Family, "^[^_]*__", ""),
      Genus_clean = str_replace(Genus, "^[^_]*__", ""),
      Species_clean = str_replace(Species, "^[^_]*__", ""),

      # 创建标签 - 优先使用Class_clean
      tax_label = case_when(
        !is.na(Class_clean) & Class_clean != "" ~ Class_clean,
        !is.na(Phylum_clean) & Phylum_clean != "" ~ paste0("p_", Phylum_clean),
        TRUE ~ "Unclassified"
      ),

      # 对重复的class名称添加序号
      tax_label = ifelse(duplicated(tax_label) | duplicated(tax_label, fromLast = TRUE),
                         paste0(tax_label, "_", ave(seq_along(tax_label), tax_label, FUN = seq_along)),
                         tax_label),

      # 短标签
      short_label = ifelse(nchar(tax_label) > 20,
                           paste0(substr(tax_label, 1, 17), "..."),
                           tax_label),

      # 图形标签 - 只使用class名称，去掉OTU_ID
      plot_label = short_label
    )

  # 计算各组的平均相对丰度
  top_otu_ids <- df_top$OTU_ID

  otu_top <- otu[, colnames(otu) %in% top_otu_ids, drop = FALSE]

  otu_with_group <- cbind(
    sample_id = rownames(otu_top),
    group = meta$treatment[match(rownames(otu_top), meta$`sample-id`)],
    otu_top
  ) %>% as.data.frame()

  # 转换为长格式
  otu_long <- otu_with_group %>%
    pivot_longer(
      cols = all_of(top_otu_ids),
      names_to = "OTU_ID",
      values_to = "relative_abundance"
    ) %>%
    mutate(
      relative_abundance = as.numeric(relative_abundance) * 100
    )

  # 计算每组的平均相对丰度
  otu_summary <- otu_long %>%
    group_by(group, OTU_ID) %>%
    summarise(
      mean_abundance = mean(relative_abundance, na.rm = TRUE),
      sd_abundance = sd(relative_abundance, na.rm = TRUE),
      n_samples = n(),
      se_abundance = sd_abundance / sqrt(n_samples),
      .groups = "drop"
    ) %>%
    filter(group %in% groups)

  # 合并分类学信息
  otu_summary <- otu_summary %>%
    left_join(df_top %>% select(OTU_ID, plot_label, tax_label, average, 
                                contribution_prop, cumulative_prop), 
              by = "OTU_ID")

  # 按贡献值排序
  otu_order <- df_top %>% arrange(desc(average)) %>% pull(OTU_ID)
  otu_summary$OTU_ID <- factor(otu_summary$OTU_ID, levels = otu_order)
  otu_summary$plot_label <- factor(otu_summary$plot_label, 
                                   levels = unique(otu_summary$plot_label[order(otu_summary$OTU_ID)]))

  # 输出统计信息
  cat("  总OTU数:", nrow(df), "\n")
  cat("  筛选OTU数:", nrow(df_top), "\n")
  cat("  累计贡献比例:", round(max(df_top$cumulative_prop) * 100, 1), "%\n")
  cat("  贡献值范围:", round(min(df_top$average), 4), "-", round(max(df_top$average), 4), "\n")
  cat("  贡献比例范围:", round(min(df_top$contribution_prop) * 100, 2), "%-", 
      round(max(df_top$contribution_prop) * 100, 2), "%\n")

  # 显示前3个重要OTU
  cat("  前3个重要OTU:\n")
  for (i in 1:min(3, nrow(df_top))) {
    cat("    ", i, ". ", df_top$OTU_ID[i], " (", df_top$tax_label[i], "): ", 
        "贡献值=", round(df_top$average[i], 4), 
        ", 贡献比例=", round(df_top$contribution_prop[i] * 100, 2), "%\n", sep = "")
  }

  # 保存结果
  write.csv(df_top, paste0("SIMPER_CORRECTED_", comp_name, "_top_OTUs.csv"), 
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(otu_summary, paste0("SIMPER_CORRECTED_", comp_name, "_abundances.csv"), 
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(df, paste0("SIMPER_CORRECTED_", comp_name, "_full_results.csv"), 
            row.names = FALSE, fileEncoding = "UTF-8")

  return(list(
    comp_name = comp_name,
    groups = groups,
    df_top = df_top,
    otu_summary = otu_summary,
    full_results = df
  ))
}

# ======================
# 5. 运行所有比较
# ======================

all_results <- list()
summary_table <- data.frame()

for (comp_name in names(comparisons)) {
  groups <- comparisons[[comp_name]]
  result <- run_corrected_simper(comp_name, groups[1], groups[2])
  all_results[[comp_name]] <- result

  # 更新汇总表
  df_top <- result$df_top
  otu_summary <- result$otu_summary

  summary_table <- rbind(summary_table, data.frame(
    Comparison = comp_name,
    Group1 = groups[1],
    Group2 = groups[2],
    Total_OTUs = nrow(result$full_results),
    Filtered_OTUs = nrow(df_top),
    Cumulative_Contribution = round(max(df_top$cumulative_prop) * 100, 1),
    Max_Contribution_Value = round(max(df_top$average), 4),
    Min_Contribution_Value = round(min(df_top$average), 4),
    Max_Abundance = round(max(otu_summary$mean_abundance, na.rm = TRUE), 2),
    Min_Abundance = round(min(otu_summary$mean_abundance[otu_summary$mean_abundance > 0], na.rm = TRUE), 2),
    stringsAsFactors = FALSE
  ))
}

# 保存汇总表
write.csv(summary_table, "SIMPER_CORRECTED_summary_all_comparisons.csv", 
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n=== 所有比较完成 ===\n")
print(summary_table)

# ======================
# 6. 计算统一的坐标轴范围
# ======================

max_abundance <- max(summary_table$Max_Abundance, na.rm = TRUE) * 1.1
max_contribution <- max(summary_table$Max_Contribution_Value, na.rm = TRUE) * 1.1

cat("\n=== 统一的坐标轴范围 ===\n")
cat("相对丰度: 0 -", round(max_abundance, 2), "%\n")
cat("贡献值: 0 -", round(max_contribution, 4), "\n")

# ======================
# 7. 为每个比较创建图表（最终修正版）
# ======================

# 1. 计算统一且合理的右轴上限（基于你的实际数据 0-11%）
unified_right_max <- 15  # 设定为15%，覆盖你的11%最大值并提供适度空间

# 2. 安全计算主坐标轴上限：取【最大丰度】和【红点所需高度】中的最大值
# 这样既能画全柱子，也能保证红点永远不会超出画面
red_dot_max_height <- unified_right_max # 红点在主轴的映射高度
safe_max_abundance <- max(max_abundance, red_dot_max_height) 

cat("\n=== 安全绘图范围设定 ===\n")
cat("统一右轴(贡献比例)上限:", unified_right_max, "%\n")
cat("安全主坐标轴(相对丰度)上限:", round(safe_max_abundance, 1), "%\n")


# 然后正常绘图...
# 3. 计算缩放因子：主轴每1%对应右轴多少%
scale_factor_final <- safe_max_abundance / unified_right_max

for (comp_name in names(all_results)) {
  result <- all_results[[comp_name]]
  df_top <- result$df_top
  otu_summary <- result$otu_summary
  groups <- result$groups
  
  # 在第7部分的循环内，绘图前加入：
  otu_summary$group <- factor(otu_summary$group, levels = groups)
  # 准备点图数据
  point_data <- df_top %>%
    mutate(
      OTU_ID = factor(OTU_ID, levels = levels(otu_summary$OTU_ID)),
      plot_label = factor(plot_label, levels = levels(otu_summary$plot_label))
    )
  
  # --- 创建图表 ---
  p <- ggplot() +
    # A. 柱状图 (相对丰度)
    geom_bar(data = otu_summary,
             aes(x = plot_label, y = mean_abundance, fill = group),
             stat = "identity", 
             position = position_dodge(width = 0.7), 
             width = 0.6,
             alpha = 0.85) +
    
    # B. 红点 (贡献比例) - 核心公式
    # y = 贡献比例(%) * 缩放因子
    geom_point(data = point_data,
               aes(x = plot_label, 
                   y = contribution_prop * 100 * scale_factor_final),
               size = 3,           # 正常点大小
               color = "#BD1100",    # 红色
               shape = 16,           # 实心圆
               alpha = 0.9,
               show.legend = FALSE) +
    
    # C. 双Y轴设置
    scale_y_continuous(
      name = "Relative Abundance (%)",       # 左轴名称
      limits = c(0, safe_max_abundance),     # 安全上限
      expand = expansion(mult = c(0, 0.05)), # 顶部留5%空白
      # 右轴定义：主轴值 / 缩放因子 = 右轴值
      sec.axis = sec_axis(~./scale_factor_final, 
                          name = "Contribution (%)", # 右轴名称
                          breaks = scales::pretty_breaks(n = 5))
    ) +
    
    # D. 颜色与标签
    # ⭐ 锁死图例顺序
   
    scale_fill_manual(
      values = group_colors, 
      name = NULL,                    # 去掉标题
      breaks = groups,          # 🚨 关键：强制图例按此顺序排列
      drop = FALSE                    # 保险：即使某组缺失也保留位置
    ) +
    # ⭐ 修改点2：去掉图表标题
    labs(title = NULL, x = NULL) +
    
    # E. 主题美化 (恢复正常合理的大小)
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 16),
      axis.title.y.right = element_text(size = 16),
      axis.title.x = element_blank(), # 移除X轴标题
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # --- 保存图像 (修正离谱的尺寸参数) ---
  # 宽度随OTU数量微调，高度固定，单位是英寸(Inch)
  plot_width <- max(8, nrow(df_top) * 0.5) 
  plot_height <- 6 
  
  ggsave(paste0("SIMPER_FINAL_", comp_name, "_plot.png"), 
         p, 
         width = plot_width, 
         height = plot_height, 
         dpi = 400, 
         limitsize = FALSE)
  
  print(p)
}

# ======================
# 8. 创建关键指标解释
# ======================

cat("\n=== SIMPER结果关键指标解释 ===\n")
cat("1. average: 贡献值，表示该OTU对组间差异的实际贡献大小\n")
cat("2. contribution_prop: 贡献比例，= average / 总average，表示该OTU的相对重要性\n")
cat("3. cumulative_prop: 累计贡献比例，表示前N个OTU累计的贡献比例\n")
cat("\n筛选标准: contribution_prop > 0.02 (贡献比例 > 2%)\n")
cat("图表说明:\n")
cat("  - 柱状图: 各组的平均相对丰度(%)\n")
cat("  - 点图: 贡献值(右Y轴)，点大小=贡献比例\n")
cat("  - 点颜色: 统一为#BD1100\n")
cat("  - 柱状图颜色: Control=#82C262, PM=#5BA2CC, OA=#F4D166, PMA=#F17060, Fe=#AC6A9F\n")
