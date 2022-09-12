# Analysis of BTX623 6/30/22 protoplasts GPS2 and GPS1-NR
library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(agricolae)
library(emmeans)
library(ggbeeswarm)
library(stringr)
library(stats)
library(knitr)


A1_before_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/A1_before_Results.csv")
A2_before_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/A2_before_Results.csv")
A3_before_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/A3_before_Results.csv")
A4_before_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/A4_before_Results.csv")
B1_before_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/B1_before_Results.csv")
B2_before_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/B2_before_Results.csv")
B3_before_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/B3_before_Results.csv")
B4_before_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/B4_before_Results.csv")
A1_after_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/A1_after_Results.csv")
A2_after_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/A2_after_Results.csv")
A3_after_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/A3_after_Results.csv")
A4_after_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/A4_after_Results.csv")
B1_after_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/B1_after_Results.csv")
B2_after_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/B2_after_Results.csv")
B3_after_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/B3_after_Results.csv")
B4_after_Results <- read_csv("/Volumes/Seagate/2022_6_30_BTX623_GPS2_GPS1NR/data/B4_after_Results.csv")

BTX_6_30_22 <- rbind(A1_before_Results, A1_after_Results, A2_before_Results, A2_after_Results, A3_before_Results, A3_after_Results, A4_before_Results, A4_after_Results, B1_before_Results, B1_after_Results, B2_before_Results, B1_after_Results, B2_before_Results, B2_after_Results, B3_before_Results, B3_after_Results, B4_before_Results, B4_after_Results)
BTX_6_30_22 <- subset(BTX_6_30_22, Area <= 300 & Area > 0)

BTX_6_30_22$TIME <- "value"
BTX_6_30_22 <- BTX_6_30_22 %>% mutate(TIME = case_when(
  grepl(pattern = "Result of Before", x = Label) ~ "-10 min",
  grepl(pattern = "slice0001_Result of After", x = Label) ~ "0h",
  grepl(pattern = "slice0002_Result of After", x = Label) ~ "2h",
  grepl(pattern = "slice0003_Result of After", x = Label) ~ "4h",
  grepl(pattern = "slice0004_Result of After", x = Label) ~ "6h",
  grepl(pattern = "slice0005_Result of After", x = Label) ~ "8h",
  grepl(pattern = "slice0006_Result of After", x = Label) ~ "10h",
  grepl(pattern = "slice0007_Result of After", x = Label) ~ "12h",
  grepl(pattern = "slice0008_Result of After", x = Label) ~ "14h",
  grepl(pattern = "slice0009_Result of After", x = Label) ~ "16h",
  grepl(pattern = "slice0010_Result of After", x = Label) ~ "18h",
  grepl(pattern = "slice0011_Result of After", x = Label) ~ "20h",
  grepl(pattern = "slice0012_Result of After", x = Label) ~ "22h",
  grepl(pattern = "slice0013_Result of After", x = Label) ~ "24h",
))

BTX_6_30_22$TREATMENT <- "value"  
BTX_6_30_22 <- BTX_6_30_22 %>% mutate(TREATMENT = case_when(
  grepl(pattern = "GA3_A1", x = Label) ~ "0uM GA3",
  grepl(pattern = "GA3_A2", x = Label) ~ "1uM GA3",
  grepl(pattern = "GA3_A3", x = Label) ~ "10uM GA3",
  grepl(pattern = "GA3_A4", x = Label) ~ "100uM GA3",
  grepl(pattern = "GA3_B1", x = Label) ~ "0uM GA3",
  grepl(pattern = "GA3_B2", x = Label) ~ "1uM GA3",
  grepl(pattern = "GA3_B3", x = Label) ~ "10uM GA3",
  grepl(pattern = "GA3_B4", x = Label) ~ "100uM GA3",
))

BTX_6_30_22$Construct <- "value"  
BTX_6_30_22 <- BTX_6_30_22 %>% mutate(Construct = case_when(
  grepl(pattern = "GA3_A", x = Label) ~ "GPS2",
  grepl(pattern = "GA3_B", x = Label) ~ "GPS1-NR",
))

colnames(BTX_6_30_22)[4] <- "Ratio"
BTX_6_30_22$Ratio <- as.numeric(BTX_6_30_22$Ratio)

data_summary <- BTX_6_30_22 %>% group_by(TIME, TREATMENT, Construct) %>% summarise(n = n(), mean = mean(Ratio), sd = sd(Ratio)) %>% mutate(sem = sd/sqrt(n), CI_lower = mean + qt((1-0.95)/2, n-1) * sem, CI_upper = mean - qt((1-0.95)/2, n-1) * sem) 

data_summary <- within(data_summary, {
  TIME <- factor(TIME, levels = c("-10 min", "0h", "2h", "4h", "6h", "8h", "10h", "12h", "14h", "16h", "18h", "20h", "22h", "24h"))
  TREATMENT <- factor(TREATMENT, levels = c("0uM GA3", "1uM GA3", "10uM GA3", "100uM GA3"))})
  
plottitle <- "Sorghum BTX623 GPS2 and GPS1-NR FRET Response to GA3"
ggplot(data_summary, aes(x = TIME, y = mean, color = Construct, group = TREATMENT, shape = TREATMENT)) + geom_errorbar(aes(ymin = mean-sem, ymax = mean+sem), width = 0.1, size = 1, position = position_dodge(width = 0.4)) + 
  geom_point(lwd = 4, position = position_dodge(width = 0.4)) + labs(x = "TIME(minutes)", y = "Ratio (FRET/CFP)", title = str_wrap(plottitle, 33)) +  theme(axis.text = element_text(size = 22)) + theme(axis.title = element_text(size = 28)) + theme(plot.title = element_text(size = 36)) + theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 22)) + scale_color_manual(values = c("GPS1-NR" = "orange", "GPS2" = "blue")) + theme(plot.title = element_text(hjust = 0.25))

# Do statistical analysis for GPS2 construct first
GPS2 <- subset(BTX_6_30_22, Construct == "GPS2")
BTX_GPS2_2way_aov <- aov(Ratio ~ TREATMENT + TIME, data = GPS2)
BTX_GPS2_2way_aov_treatment_x_time <- aov(Ratio ~ TREATMENT*TIME, data = GPS2)
summary(BTX_GPS2_2way_aov)
summary(BTX_GPS2_2way_aov_treatment_x_time)


# check the gets model fit using AIC
AIC(BTX_GPS2_2way_aov, BTX_GPS2_2way_aov_treatment_x_time)
# says the better fit is the standard 2way anova, not the one with time as an interaction variable, but sam wants the interaction variable pval included
# check homoscedasticity or assumption of equal or similar variances in different groups being compared
plot(BTX_GPS2_2way_aov)
# homoscedasticity looks okay I think
# Now lets do the post hoc test
TukeyHSD(BTX_GPS2_2way_aov)
# Try plotting this in a visueally pleasing way
tukey.plot.aov <- aov(Ratio~TREATMENT:TIME, data = GPS2)
tukey.plot.test <- TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)
# not really worth it, instead lets just print the anova table
summary.aov(BTX_GPS2_2way_aov)
kable(BTX_GPS2_2way_aov, digits = 2)
#post hoc GPS2
HSD.test(BTX_GPS2_2way_aov, "TREATMENT", group = T, console = T)

# Do statistical analysis for GPS1-NR construct 
GPS1_NR <- subset(BTX_6_30_22, Construct == "GPS1-NR")
BTX_GPS1NR_2way_aov <- aov(Ratio ~ TREATMENT+TIME, data = GPS1_NR)
summary(BTX_GPS1NR_2way_aov)
HSD.test(BTX_GPS1NR_2way_aov, "TREATMENT", group = T, console = T)

# Now plot the number of protoplast nuclei
nplottitle <- "Sorghum BTX623 protoplast nuclei numbers"
ggplot(data_summary, aes(x = TIME, y = n, color = Construct, group = TREATMENT, shape = TREATMENT)) + geom_point() + labs(x = "Time", y = "N protoplasts", title = str_wrap(nplottitle, 30)) + scale_color_manual(values = c("GPS1-NR" = "darkorange2", "GPS2" = "blue")) +
  theme(axis.text = element_text(size = 12)) + theme(axis.title = element_text(size = 26)) + theme(plot.title = element_text(size = 32)) + theme(legend.text = element_text(size = 18)) + theme(legend.title = element_text(size = 18)) + theme(plot.title = element_text(hjust = 0.25)) + theme(panel.background = element_rect(fill = "papayawhip")) + scale_y_continuous(expand = c(0, 0), limits = c(0, 90))





