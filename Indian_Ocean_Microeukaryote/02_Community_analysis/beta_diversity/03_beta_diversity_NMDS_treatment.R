#Bray-Curtis distance matrix to NMDS plot
library(vegan)
library(readr)
library(dplyr)
library(ggplot2)

# ======================
# 1. 读入距离矩阵
# ======================
dist <- read_tsv("C:/Users/26355/Desktop/xxx/beta/NMDS/distance-matrix.tsv")

dist_mat <- as.matrix(dist[,-1])
rownames(dist_mat) <- dist[[1]]

dist_obj <- as.dist(dist_mat)

# ======================
# 2. NMDS分析
# ======================
nmds <- metaMDS(dist_obj, k = 2)

# 提取坐标
points <- as.data.frame(nmds$points)
points$SampleID <- rownames(points)

# ======================
# 3. 合并metadata
# ======================
meta <- read_tsv("C:/Users/26355/Desktop/xxx/beta/NMDS/metadata.tsv") %>%
  rename(SampleID = `sample-id`)

data <- left_join(points, meta, by = "SampleID")

# 4. 设置分组
# ======================
data$treatment <- factor(data$treatment,
                         levels = c("Control", "PM", "A", "PMA", "Fe"))

# 计算stress值
stress_value <- nmds$stress
stress_label <- paste0("Stress: ", round(stress_value, 4))

# ======================
# 5. 画图（应用颜色格式）
# ======================
p <- ggplot(data, aes(MDS1, MDS2, color = treatment)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(color = treatment),
               level = 0.995,
               linetype = 2,
               linewidth = 0.8) +
  scale_color_manual(values = c(
    "Control" = "#00AA5A",
    "PM" = "#00A6CA",
    "A" = "#EEB53C",
    "PMA" = "#C52C2A",
    "Fe" = "#9F7FE5"
  )) +
  # 添加stress值标注
  annotate("text", 
           x = -Inf, 
           y = -Inf,
           label = stress_label,
           size = 3.5,  # 字体稍小
           hjust = -0.1,   # 左对齐
           vjust = -1,   # 底部对齐
           color = "black") +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    legend.position = "right"
  ) +
  labs(
    x = "MDS1",
    y = "MDS2",
    color = "Treatment"
  )

print(p)

# ======================
# 6. 保存图片
# ======================
# 保存到当前工作目录
ggsave("NMDS_treatment_stress.png", plot = p,
       width = 6, height = 5, dpi = 300)

# 保存到指定路径
ggsave("C:/Users/26355/Desktop/xxx/beta/NMDS/NMDS_treatment_stress.png",
       plot = p,
       width = 6, height = 5, dpi = 300)