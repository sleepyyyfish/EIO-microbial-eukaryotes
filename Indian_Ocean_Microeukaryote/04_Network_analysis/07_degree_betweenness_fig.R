# Figure5b_完整分析.R

# 描述：计算和绘制标准化度与介数中心性的双对数关系图

setwd("C:/Users/26355/Desktop/xxx/network analysis")
# 加载包
library(igraph)
library(tidyverse)
library(ggplot2)
library(patchwork)

# 构建网络
read_network <- function(node_file, edge_file){
  
  nodes <- read.csv(node_file, stringsAsFactors = FALSE)
  edges <- read.csv(edge_file, stringsAsFactors = FALSE)
  
  colnames(nodes)[1] <- "id"
  
  # ⭐ 关键：过滤负相关
  edges <- edges %>%
    filter(weight > 0)
  
  g <- graph_from_data_frame(d = edges,
                             vertices = nodes,
                             directed = FALSE)
  
  # 删除孤立点
  g <- delete.vertices(g, degree(g) == 0)
  
  return(g)
}

#读取五个网络
g_control <- read_network("Control_nodes.csv", "Control_edges.csv")
g_PM      <- read_network("PM_nodes.csv", "PM_edges.csv")
g_A      <- read_network("A_nodes.csv", "A_edges.csv")
g_PMA     <- read_network("PMA_nodes.csv", "PMA_edges.csv")
g_Fe      <- read_network("Fe_nodes.csv", "Fe_edges.csv")

# 步骤1：计算函数
calc_degree_between <- function(g){
  N <- vcount(g)
  deg <- degree(g)
  norm_deg <- deg / (N - 1)
  btw <- betweenness(g, normalized = FALSE)
  
  df <- data.frame(
    node = V(g)$name,
    degree = deg,
    norm_degree = norm_deg,
    betweenness = btw
  )
  
  df <- df %>%
    filter(norm_degree > 0, betweenness > 0) %>%
    mutate(
      log_deg = log10(norm_degree),
      log_btw = log10(betweenness)
    )
  
  return(df)
}

# 步骤2：计算所有组
db_control <- calc_degree_between(g_control) %>% mutate(group = "Control")
db_PM <- calc_degree_between(g_PM) %>% mutate(group = "PM")
db_A <- calc_degree_between(g_A) %>% mutate(group = "A")
db_PMA <- calc_degree_between(g_PMA) %>% mutate(group = "PMA")
db_Fe <- calc_degree_between(g_Fe) %>% mutate(group = "Fe")

db_all <- rbind(db_control, db_PM, db_A, db_PMA, db_Fe)
db_all$group <- factor(db_all$group, levels = c("Control", "PM", "A", "PMA", "Fe"))

# 步骤3：计算统计
cor_stats <- db_all %>%
  group_by(group) %>%
  summarise(
    n = n(),
    cor = cor(log_deg, log_btw, method = "pearson"),
    p_value = cor.test(log_deg, log_btw, method = "pearson")$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "NS"
    )
  )

# 步骤4：绘图 - 只显示r值
# 定义自定义颜色
group_colors <- c(
  "Control" = "#00AA5A",  # 绿色
  "PM" = "#00A6CA",       # 蓝色
  "A" = "#EEB53C",        # 橙色
  "PMA" = "#C52C2A",      # 红色
  "Fe" = "#9F7FE5"        # 紫色
)

p <- ggplot(db_all, aes(x = log_deg, y = log_btw, color = group)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  facet_wrap(~group) +
  scale_color_manual(values = group_colors) +
  # 只显示r值
  geom_text(
    data = cor_stats,
    aes(x = Inf, y = -Inf,  # x=最大值，y=最小值
        label = sprintf("r = %.3f", cor),
        hjust = 1.1,  # 水平：1=右对齐，>1=更靠右
        vjust = -0.5), # 垂直：0=底对齐，<0=更靠下
    color = "black", size = 5, inherit.aes = FALSE
  ) +
  facet_wrap(~group, nrow = 1) +
  theme_bw(base_size = 16) +
  labs(
    x = expression(log[10]*"(Normalized degree)"),
    y = expression(log[10]*"(Betweenness centrality)"),
  ) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold", size = 16),
    axis.title.x = element_text(size = 16),  # x轴标题大小
    axis.title.y = element_text(size = 16)  # y轴标题大小
  )

print(p)

# 显示结果

print(cor_stats)

# 保存
ggsave("Figure5b_degree_betweenness.png", width = 20, height = 5, dpi = 300)
cat("分析完成！图形已保存为Figure5b_degree_betweenness.png\n")
cat("\n相关性统计结果：\n")
print(cor_stats)