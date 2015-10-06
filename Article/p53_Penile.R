library(foreign)
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
Data <- read.dta("p53_tma_long.dta")

# Transforming the Data
Data <- Data %>%
  mutate(
    p53_z = log10(p53 + 1),
    ijp53_z = log10(ijp53 + 1)
  )

# Descriptive Statistics
No_Cases <- length(unique(Data$caseid))
No_Spots <- nrow(Data)
p53_Mean <- round(mean(Data$p53, na.rm = TRUE), 1)
# p53 visual estimation
p53_SD <- round(sd(Data$p53, na.rm = TRUE), 1)
p53_Min <- round(min(Data$p53, na.rm = TRUE), 1)
p53_Max <- round(max(Data$p53, na.rm = TRUE), 1)
p53_Mean <- round(mean(Data$p53, na.rm = TRUE), 1)
p53_SD <- round(sd(Data$p53, na.rm = TRUE), 1)
p53_Min <- round(min(Data$p53, na.rm = TRUE), 1)
p53_Max <- round(max(Data$p53, na.rm = TRUE), 1)
# p53 digital estimation
p53ij_Mean <- round(mean(Data$ijp53, na.rm = TRUE), 1)
p53ij_SD <- round(sd(Data$ijp53, na.rm = TRUE), 1)
p53ij_Min <- round(min(Data$ijp53, na.rm = TRUE), 1)
p53ij_Max <- round(max(Data$ijp53, na.rm = TRUE), 1)

# Inferential Statistics
# p53 visual vs digital estimation
p53_Wilcox <- format(wilcox.test(Data$p53, Data$ijp53)$p.value, digits = 2)
p53_Spear <- cor.test(Data$p53, Data$ijp53, method = "spearman")
p53_Spear_Est <- round(p53_Spear$estimate, 2)
p53_Spear_P <- format(p53_Spear$p.value, digits = 2)

# Reshaping for counting subtype and grade by cases
Data_Wide <- plyr::ddply(Data, .(caseid, subtype), summarize, grade = max(grade, na.rm = TRUE))

# Reshaping data for Figure 2
Data_Long <- Data %>%
  select(p53_z, ijp53_z) %>%
  gather(method, value, p53_z:ijp53_z)

# Data for Table 1
## Histologic subtype
T1_s <- Data_Wide %>%
  group_by(subtype) %>%
  summarize(N = n()) %>%
  arrange(desc(N))
## Histologic grade
T1_hg <- Data_Wide %>%
  group_by(grade) %>%
  summarize(N = n())
## Visual & digital estimation by histologic subtype
T1_p53_s <- Data %>%
  group_by(subtype) %>%
  summarize(
    N = n(),
    Mean_V = round(mean(p53, na.rm = TRUE), 1),
    SD_V = round(sd(p53, na.rm = TRUE), 1),
    Mean_D = round(mean(ijp53, na.rm = TRUE), 1),
    SD_D = round(sd(ijp53, na.rm = TRUE), 1)
  ) %>%
  arrange(desc(N))
## Visual & digital estimation by histologic grade
T1_p53_hg <- Data %>%
  group_by(grade) %>%
  filter(!is.na(grade)) %>%
  summarize(
    N = n(),
    Mean_V = round(mean(p53, na.rm = TRUE), 1),
    SD_V = round(sd(p53, na.rm = TRUE), 1),
    Mean_D = round(mean(ijp53, na.rm = TRUE), 1),
    SD_D = round(sd(ijp53, na.rm = TRUE), 1)
  )
## Krukal-Wallis tests
KW_p53_s <- format(kruskal.test(p53 ~ subtype, Data)$p.value, digits = 2)
KW_ijp53_s <- format(kruskal.test(ijp53 ~ subtype, Data)$p.value, digits = 2)
KW_p53_hg <- format(kruskal.test(p53 ~ grade, Data)$p.value, digits = 2)
KW_ijp53_hg <- format(kruskal.test(ijp53 ~ grade, Data)$p.value, digits = 2)

## Figure 2
library(cowplot)
plot1 <- ggplot(Data_Long, aes(method, value, fill = method)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  guides(fill = FALSE) +
  scale_x_discrete(labels = c("Visual Estimation", "Digital Estimation")) +
  ylab("Positive Cells, log10(%)") + xlab("") +
  theme_cowplot(font_size = 10)
plot2 <- ggplot(Data, aes(p53_z, ijp53_z)) +
  geom_point(shape = 1) +
  geom_smooth(method = lm, se = FALSE) +
  ylab("Visual Estimation, log10(%)") +
  xlab("Digital Estimation, log10(%)") +
  theme_cowplot(font_size = 10)
plot_grid(plot1, plot2, labels = c("A", "B"))
ggsave("Fig2.png", height = 3)
