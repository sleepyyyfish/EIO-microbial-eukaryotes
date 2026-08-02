# ======================
# 0. 加载包（缺哪个就先install.packages）
# ======================
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# ======================
# 1. 读入距离矩阵
# ======================
dist <- read_tsv("C:/Users/26355/Desktop/xxx/beta/平均相异度/distance-matrix.tsv")

dist_mat <- as.matrix(dist[,-1])
rownames(dist_mat) <- dist[[1]]

dist_df <- as.data.frame(dist_mat)

# ======================
# 2. 转成长表
# ======================
dist_long <- dist_df %>%
  rownames_to_column("Sample1") %>%
  pivot_longer(-Sample1, names_to = "Sample2", values_to = "distance") %>%
  filter(Sample1 < Sample2)

# ======================
# 3. 读metadata
# ======================
meta <- read_tsv("C:/Users/26355/Desktop/毕设111/aaa重来一遍/beta/平均相异度/metadata.tsv") %>%
  rename(SampleID = `sample-id`) %>%
  mutate(
    treatment = factor(treatment, 
                       levels = c("Control", "PM", "A", "PMA", "Fe")),
    day = paste0("d", day)
  ) %>%
  mutate(day = factor(day, 
                      levels = c("d1", "d3", "d5", "d9", "d16", "d20")))

# 定义组比较顺序
treatment_levels <- c("Control", "PM", "A", "PMA", "Fe")
pair_levels <- apply(combn(treatment_levels, 2), 2, 
                     function(x) paste(x[1], "vs", x[2]))

# ======================
# 4. 合并信息
# ======================
dist_long <- dist_long %>%
  left_join(meta, by = c("Sample1" = "SampleID")) %>%
  rename(treatment1 = treatment, day1 = day) %>%
  left_join(meta, by = c("Sample2" = "SampleID")) %>%
  rename(treatment2 = treatment, day2 = day)

# ======================
# 5. 统一pair（关键修复部分）
# ======================
dist_long <- dist_long %>%
  filter(treatment1 != treatment2) %>%
  mutate(
    # 将treatment转换为数值索引
    idx1 = as.numeric(treatment1),
    idx2 = as.numeric(treatment2),
    # 按照预定义的treatment_levels顺序排序
    g1 = ifelse(idx1 <= idx2, as.character(treatment1), as.character(treatment2)),
    g2 = ifelse(idx1 <= idx2, as.character(treatment2), as.character(treatment1)),
    pair = paste(g1, "vs", g2)
  ) %>%
  mutate(
    pair = factor(pair, levels = pair_levels)
  ) %>%
  select(-idx1, -idx2)  # 移除临时列

# ======================
# 6. (d) 组间平均距离
# ======================
dist_group <- dist_long %>%
  group_by(g1, g2) %>%
  summarise(mean_dist = mean(distance), .groups = "drop") %>%
  mutate(
    g1 = factor(g1, levels = treatment_levels),
    g2 = factor(g2, levels = treatment_levels)
  )

p_d <- ggplot(dist_group, aes(x = g1, y = g2, fill = mean_dist)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_dist, 2)), size = 4) +
  scale_fill_gradient(low = "#F5F5F5", high = "#3C5488") +
  theme_classic() +
  labs(x = "", y = "", fill = "Distance")

print(p_d)

ggsave("C:/Users/26355/Desktop/毕设111/aaa重来一遍/beta/平均相异度/d_heatmap.png",
       plot = p_d, width = 5, height = 4, dpi = 300)

# ======================
# 7. (e) 按时间分
# ======================
dist_day_group <- dist_long %>%
  filter(day1 == day2) %>%
  group_by(day = day1, pair) %>%
  summarise(mean_dist = mean(distance), .groups = "drop") %>%
  mutate(
    day = factor(day, levels = c("d1", "d3", "d5", "d9", "d16", "d20")),
    pair = factor(pair, levels = pair_levels)
  )

p_e <- ggplot(dist_day_group, aes(x = pair, y = day, fill = mean_dist)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(mean_dist, 2)), size = 3) +
  scale_fill_gradient(low = "#F5F5F5", high = "#3C5488") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(x = "Group comparison", y = "Day")

print(p_e)

ggsave("C:/Users/26355/Desktop/xxx/beta/平均相异度/e_heatmap.png",
       plot = p_e, width = 7, height = 4, dpi = 300)

# remember to add ** using permanova results by yourself