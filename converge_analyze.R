library(readxl)
library(ggplot2)
library(mgcv)

# 1. 读入原始数据
# 请确保路径指向你的实际文件
df <- read_excel("C:/Paper and funding/2024 Autotrophy niche/Converge分析.xlsx")

# 2. 重命名列
colnames(df) <- c("SAM_depth", "Equal_PD_AP_AOA", "PNM_depth", "Equal_PD_NP_NOB")

# 3. 添加时间列 (1 到 1000)
df$time <- 1:nrow(df)

# 4. 计算两个深度差值
df$delta_SAM_PD <- df$SAM_depth - df$Equal_PD_AP_AOA
df$delta_PNM_PD <- df$PNM_depth - df$Equal_PD_NP_NOB

# 5. 合并数据：转换成长格式
df_long <- rbind(
  data.frame(time = df$time, delta = df$delta_SAM_PD, group = "SAM – Equal PD"),
  data.frame(time = df$time, delta = df$delta_PNM_PD, group = "PNM – Equal PD")
)

# 控制图例顺序
df_long$group <- factor(df_long$group, levels = c("SAM – Equal PD", "PNM – Equal PD"))

# 设置颜色映射
color_map <- c("SAM – Equal PD" = "#4B65AF", "PNM – Equal PD" = "#AF4B4B")

# 6. 绘图
ggplot(df_long, aes(x = time, y = delta, color = group, fill = group)) +
  # 散点图
  geom_point(size = 0.8, alpha = 0.5, color = "gray60") +
  # 全局 GAM 拟合曲线 (基于所有 1000 天数据)
  geom_smooth(method = "gam", formula = y ~ s(x), linewidth = 1.2, alpha = 0.3) +
  # 颜色设置
  scale_color_manual(values = color_map) +
  scale_fill_manual(values = color_map) +
  # 坐标轴标签逻辑修改
  scale_x_continuous(
    # 每 30 天一个刻度，从第 1 天开始：1, 31, 61...
    breaks = seq(1, 301, by = 30), 
    # 将刻度值转化为标签：(1-1)/30=0, (31-1)/30=1...
    labels = function(x) (x - 1) / 30,
    expand = c(0, 0)
  ) +
  # 关键：使用 coord_cartesian 进行“放大”，不影响 geom_smooth 的全局计算
  coord_cartesian(xlim = c(1, 301)) + 
  labs(
    title = "Convergence of SAM and PNM toward Ecological Interface",
    x = "Months after mixing", 
    y = expression(bold(Delta*"Depth (m)")),
    color = NULL, fill = NULL
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Arial", face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10, face = "bold"),
    legend.position = "top",
    aspect.ratio = 0.4
  )