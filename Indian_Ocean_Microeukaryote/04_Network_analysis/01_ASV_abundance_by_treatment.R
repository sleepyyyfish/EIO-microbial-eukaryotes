# ======================
# 0. 设置路径 + 加载包
# ======================
setwd("C:/Users/26355/Desktop/xxx/network analysis/")

library(tidyverse)

# ======================
# 1. 读入数据
# ======================

# 丰度表（行为ASV，列为样本）
otu <- read_csv("feature-table.csv")

# metadata
meta <- read_tsv("metadata.tsv")

# ======================
# 2. 整理列名
# ======================

# 确保第一列是ASV ID
colnames(otu)[1] <- "ASV"

# ======================
# 3. 宽表转长表
# ======================

otu_long <- otu %>%
  pivot_longer(
    cols = -ASV,
    names_to = "sample-id",
    values_to = "abundance"
  )

# ======================
# 4. 合并 metadata
# ======================

otu_meta <- otu_long %>%
  left_join(meta, by = "sample-id")

# ======================
# 5. 按 treatment 汇总（核心步骤！）
# ======================

otu_sum <- otu_meta %>%
  group_by(ASV, treatment) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    .groups = "drop"
  )

# ======================
# 6. 再转回宽表（每个处理一列）
# ======================

otu_wide <- otu_sum %>%
  pivot_wider(
    names_from = treatment,
    values_from = total_abundance,
    values_fill = 0
  )

# ======================
# 7. 输出结果
# ======================

write_csv(otu_wide, "ASV_abundance_by_treatment.csv")