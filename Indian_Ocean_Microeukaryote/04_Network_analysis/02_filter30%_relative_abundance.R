# ===============================
# Step 0｜路径
# ===============================
setwd("C:/Users/26355/Desktop/xxx/network analysis/")

library(tidyverse)

# ===============================
# Step 1｜读入原始feature-table
# ===============================
asv_raw <- read.table("feature-table.tsv",
                      header = TRUE,
                      row.names = 1,
                      sep = "\t",
                      check.names = FALSE,
                      comment.char = "")

# ❗检查一下
dim(asv_raw)
head(asv_raw[,1:5])

# 👉 此时应该是：
# 行 = ASV（长ID）
# 列 = 样本（M01, M02...）

# ===============================
# Step 2｜转置（变成 行=样本）
# ===============================
asv_t <- t(asv_raw)

# 检查
dim(asv_t)
rownames(asv_t)[1:5]   # 应该是 M01 M02...
colnames(asv_t)[1:5]   # 应该是 ASV长ID

# ===============================
# Step 3｜转相对丰度
# ===============================
asv_rel <- asv_t / rowSums(asv_t)

# ===============================
# Step 4｜筛选 ≥30%出现（不多不少，节点几十到几百）
# ===============================
asv_occurrence <- colSums(asv_rel > 0)
threshold <- nrow(asv_rel) * 0.3

asv_filtered <- asv_rel[, asv_occurrence >= threshold]

# ===============================
# Step 5｜保存（正确版本）
# ===============================
write.csv(asv_filtered, "ASV_filtered_rel_abundance_CORRECT.csv")