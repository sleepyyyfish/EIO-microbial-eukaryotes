# ===============================
# Step 0｜路径 + 加载包
# ===============================
setwd("C:/Users/26355/Desktop/xxx/network analysis")

library(tidyverse)
library(Hmisc)

# ===============================
# Step 1｜读入ASV数据（正确版本）
# ===============================
asv <- read.csv("ASV_filtered_rel_abundance_30.csv",
                row.names = 1,
                check.names = FALSE)

# 检查
head(rownames(asv))   # 应该是 M01 M02...
head(colnames(asv))   # 应该是 ASV长ID

# ===============================
# Step 2｜读入metadata
# ===============================
meta <- read.table("metadata.tsv",
                   header = TRUE,
                   sep = "\t",
                   check.names = FALSE)

# 重命名（关键！）
meta <- meta %>%
  rename(sample = `sample-id`)

# ===============================
# Step 3｜合并数据
# ===============================
asv$sample <- rownames(asv)

# 合并
data <- left_join(meta, asv, by = "sample")

# 修复列名问题（关键！）
colnames(data) <- make.names(colnames(data), unique = TRUE)

# 再检查
any(colnames(data) == "")
any(is.na(colnames(data)))

# 检查是否匹配成功
sum(is.na(data$treatment))  # 应该是0

# ===============================
# Step 4｜建网络函数
# ===============================
build_network <- function(df, group_name){
  
  df_mat <- df %>% 
    select(-sample, -treatment, -description, -group, -replicate, -day) %>%
    as.matrix()
  
  # 计算相关性
  cor_res <- rcorr(df_mat, type = "spearman")
  
  cor_matrix <- cor_res$r
  p_matrix <- cor_res$P
  
  # 阈值
  cor_threshold <- 0.6
  p_threshold <- 0.05
  
  edges <- which(abs(cor_matrix) > cor_threshold & p_matrix < p_threshold, arr.ind = TRUE)
  
  edge_list <- data.frame(
    source = rownames(cor_matrix)[edges[,1]],
    target   = colnames(cor_matrix)[edges[,2]],
    weight = cor_matrix[edges],
    pvalue = p_matrix[edges]
  )
  
  # 去重 - 修改这里：将from/to改为source/target
  edge_list <- edge_list[edge_list$source != edge_list$target, ]
  edge_list <- edge_list[!duplicated(t(apply(edge_list[,1:2], 1, sort))), ]
  
  nodes <- data.frame(id = colnames(df_mat))
  
  # 新增：去掉"X"前缀
  edge_list$source <- sub("^X", "", edge_list$source)
  edge_list$target <- sub("^X", "", edge_list$target)
  nodes$id <- sub("^X", "", nodes$id)
  
  # 保存
  write.csv(edge_list, paste0(group_name, "_edges.csv"), row.names = FALSE)
  write.csv(nodes, paste0(group_name, "_nodes.csv"), row.names = FALSE)
}

# ===============================
# Step 5｜按treatment分组建网络
# ===============================
groups <- unique(data$treatment)

for(g in groups){
  
  df_sub <- data %>% filter(treatment == g)
  
  build_network(df_sub, g)
}