setwd("C:/Users/26355/Desktop/xxx/network analysis")
# =========================
# 0. packages
# =========================
library(dplyr)
library(readr)
library(readxl)

# =========================
# 1. load data (以Fe为例，其他依次修改运行)
# =========================
nodes <- read_csv("Fe_nodes.csv")
tax <- read_csv("all_asv_tax.csv")
abund <- read_excel("ASV_abundance_by_treatment.xlsx")

# =========================
# 2. clean column names
# =========================
colnames(nodes)[1] <- "ASV"
nodes$ASV <- trimws(as.character(nodes$ASV))

# =========================
# 3. add taxonomy
# =========================
nodes <- nodes %>%
  left_join(tax, by = c("ASV" = "OUT_ID"))

# 处理缺失的分类信息
nodes$Phylum[is.na(nodes$Phylum)] <- "others"

# =========================
# 4. add abundance
# =========================
nodes <- nodes %>%
  left_join(abund %>% select(ASV, Fe), by = "ASV")

# 处理缺失的丰度值
nodes$Fe[is.na(nodes$Fe)] <- 0

# =========================
# 5. 创建Gephi节点文件
# =========================
gephi_Fe_nodes <- nodes %>%
  select(
    id = ASV,           # 节点ID
    phylum = Phylum,    # 门分类
    abundance = Fe  # 丰度
  ) %>%
  # 确保没有重复的行
  distinct(id, .keep_all = TRUE)

# =========================
# 6. 保存为CSV文件
# =========================
write_csv(gephi_Fe_nodes, "Gephi_Fe_nodes.csv")

cat("✅ 文件已生成: GGephi_Fe_nodes.csv\n")
cat("包含", nrow(gephi_Fe_nodes), "个节点\n")
cat("列名: id, phylum, abundance\n")