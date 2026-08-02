setwd("C:/Users/26355/Desktop/毕设111/aaa重来一遍/网络分析network analysis")

library(igraph)
library(dplyr)
library(ggplot2)

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

# ziPi计算

calc_zipi <- function(g){
  
  # 1. 模块划分
  comm <- cluster_fast_greedy(g)
  module <- membership(comm)
  V(g)$module <- module
  
  # 2. Zi计算
  Zi <- numeric(vcount(g))
  
  for(i in 1:vcount(g)){
    v <- V(g)[i]
    mod <- module[i]
    
    nodes_in_mod <- which(module == mod)
    
    # 节点在模块内的连接数
    nei <- neighbors(g, v)
    k_i <- sum(nei %in% nodes_in_mod)
    
    # 模块内所有节点的连接数
    k_s <- sapply(nodes_in_mod, function(n){
      sum(neighbors(g, V(g)[n]) %in% nodes_in_mod)
    })
    
    if(sd(k_s) == 0){
      Zi[i] <- 0
    } else {
      Zi[i] <- (k_i - mean(k_s)) / sd(k_s)
    }
  }
  
  # 3. Pi计算
  Pi <- numeric(vcount(g))
  
  for(i in 1:vcount(g)){
    v <- V(g)[i]
    k_i <- degree(g, v)
    
    if(k_i == 0){
      Pi[i] <- 0
      next
    }
    
    mods <- unique(module)
    sum_sq <- 0
    
    for(m in mods){
      nodes_in_mod <- which(module == m)
      k_im <- sum(neighbors(g, v) %in% nodes_in_mod)
      sum_sq <- sum_sq + (k_im / k_i)^2
    }
    
    Pi[i] <- 1 - sum_sq
  }
  
  # 4. 分类
  role <- ifelse(Zi > 2.5 & Pi > 0.62, "Network hub",
                 ifelse(Zi > 2.5 & Pi <= 0.62, "Module hub",
                        ifelse(Zi <= 2.5 & Pi > 0.62, "Connector",
                               "Peripheral")))
  
  df <- data.frame(
    node = V(g)$name,
    module = module,
    Zi = Zi,
    Pi = Pi,
    role = role
  )
  
  return(df)
}

# 批量五组
zipi_control <- calc_zipi(g_control)
zipi_PM      <- calc_zipi(g_PM)
zipi_A      <- calc_zipi(g_A)
zipi_PMA     <- calc_zipi(g_PMA)
zipi_Fe      <- calc_zipi(g_Fe)



#画图
zipi_control$group <- "Control"
zipi_PM$group <- "PM"
zipi_A$group <- "A"
zipi_PMA$group <- "PMA"
zipi_Fe$group <- "Fe"

zipi_all <- rbind(zipi_control, zipi_PM, zipi_A, zipi_PMA, zipi_Fe)

# 👇 添加这行代码设置因子顺序
zipi_all$group <- factor(zipi_all$group,
                         levels = c("Control", "PM", "A", "PMA", "Fe"))

# ...前面的代码保持不变，直到画图部分...

# 读取分类信息
taxonomy <- read.csv("all_asv_tax.csv", stringsAsFactors = FALSE)
colnames(taxonomy) <- c("node", "Phylum")  # 确保列名匹配

# 合并分类信息到zipi_all
zipi_all <- merge(zipi_all, taxonomy, by = "node", all.x = TRUE)

# 创建标签列：只对非Peripheral点显示Phylum
zipi_all$label <- ifelse(zipi_all$role != "Peripheral", zipi_all$Phylum, "")

# ... 前面的数据处理代码 ...

# 设置因子顺序
zipi_all$group <- factor(zipi_all$group,
                         levels = c("Control", "PM", "A", "PMA", "Fe"))

# 确保role列是因子，并包含所有可能的角色
role_levels <- c("Network hub", "Module hub", "Connector", "Peripheral")
zipi_all$role <- factor(zipi_all$role, levels = role_levels)

# 创建颜色映射
role_colors <- c("Network hub" = "#2E5A87",
                 "Module hub" = "#E41A1C",
                 "Connector" = "#DC6826",
                 "Peripheral" = "#999999")

# 创建图表
p <- ggplot(zipi_all, aes(x = Pi, y = Zi, color = role)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_vline(xintercept = 0.62, linetype = "dashed", alpha = 0.5) +
  geom_hline(yintercept = 2.5, linetype = "dashed", alpha = 0.5) +
  geom_text_repel(aes(label = label), 
                  size = 4, 
                  color = "black",
                  box.padding = 0.8,
                  point.padding = 0.3, 
                  max.time = 2,  # 减少计算时间
                  max.iter = 10000,  # 增加迭代次数
                  force = 1,  # 增加排斥力
                  min.segment.length = 0.2,  # 增加连线的最小长度
                  max.overlaps = Inf,
                  segment.size = 0.2,
                  segment.color = "grey50",
                  show.legend = FALSE) +
  facet_wrap(~group, nrow = 1) +
  scale_color_manual(values = role_colors, drop = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 16, face = "bold"),
        # 在这里添加这两行👇
        axis.title.x = element_text(size = 16),  # x轴标题大小
        axis.title.y = element_text(size = 16)) +  # y轴标题大小
  labs(x = "Among-module connectivity (Pi)",
       y = "Within-module connectivity (Zi)",
       color = "Node role")

# 显示图表
print(p)


ggsave("Zi-Pi_topological_role_network_hubs.pdf", p, width = 20, height = 5, device = cairo_pdf)

# 也保存一个高分辨率PNG备用
ggsave("Zi-Pi_topological_role_network_hubs.png", p, width = 20, height = 5, dpi = 600, bg = "white")

# 如果你想查看非Peripheral点的详细信息
non_peripheral <- zipi_all[zipi_all$role != "Peripheral", ]
print("非Peripheral点的详细信息：")
print(non_peripheral)

# 按组和角色统计
role_counts <- table(zipi_all$group, zipi_all$role)
print("各处理组的节点角色统计：")
print(role_counts)

# 也可以保存统计结果为CSV
write.csv(as.data.frame.matrix(role_counts), "node_role_statistics.csv")