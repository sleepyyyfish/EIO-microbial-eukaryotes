setwd("C:/Users/26355/Desktop/xxx/alpha")

library(tidyverse)
#把提取后单独成两列的数据长表转宽表
#以SampleID|Shannon为例
# 1. 读取数据
alpha <- read.table("alpha-diversity-Shannon.tsv",
                    header = TRUE,
                    sep = "\t",
                    stringsAsFactors = FALSE)

metadata <- read.table("metadata.tsv",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE)

# 2. 合并 metadata + alpha diversity
colnames(alpha)[1] <- "sample_id"
colnames(metadata)[1] <- "sample_id"

df <- left_join(metadata, alpha, by = "sample_id")

# 3. 创建“列名组合”（Control_rep1这种结构）
df2 <- df %>%
  mutate(group_rep = paste0(treatment, "_rep", replicate))

# 4. 转宽表：day作为行，group_rep作为列
wide <- df2 %>%
  select(day, group_rep, Shannon) %>%
  pivot_wider(
    names_from = group_rep,
    values_from = Shannon
  )

# 5. 排序行（按时间顺序）
day_order <- c("d1", "d3", "d5", "d9", "d16", "d20")

wide$day <- factor(wide$day, levels = day_order)
wide <- wide %>% arrange(day)

# 6. （可选）如果你想只保留你说的顺序列
col_order <- c(
  "Control_rep1","Control_rep2","Control_rep3",
  "PM_rep1","PM_rep2","PM_rep3",
  "A_rep1","A_rep2","A_rep3",
  "PMA_rep1","PMA_rep2","PMA_rep3",
  "Fe_rep1","Fe_rep2","Fe_rep3"
)

wide <- wide %>%
  select(day, any_of(col_order))

# 7. 输出
write.table(wide,
            file = "shannon_wide_table.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)