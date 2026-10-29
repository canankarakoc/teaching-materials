# Practical training UL
# Bacterial growth models / Virus infectivity 
# 14.05.20 Karakoc 

set.seed(140520)
# install.packages("growthrates")
# see https://tpetzoldt.github.io/growthrates/doc/Introduction.html

library(growthrates)
library(tidyverse)

# more about tidywerse (dplyr and ggplot2)
# https://www.tidyverse.org/

# set working directory
# run from the folder this script is in

# read data
# you can play with the growth_data1 
# it is identical with one difference which is Local_line (swapped bacteria and phages or sympatric/allopatric phages)

# OD data from the plate reader 
growth_data2  <- read.table("growth_data2.csv", header=T, sep=";", dec=",")

# your plate design
plate_design2 <- read.table("plate_design2.csv", header=T, sep=";", dec=",")

merged_data <- growth_data2 %>%
  gather("Well", "OD", -Time) %>%
  full_join(plate_design2, by = "Well") %>%
  filter(SampleID != "blank") %>%
  group_by(SampleID, Time, Sample_line, Treatment, Phage_line, Local_line) %>%
  summarise(OD = mean(OD))
  
  
ggplot(merged_data, aes(x = Time, y = OD))+
  geom_line(color = "darkred")+
  facet_wrap(SampleID~., ncol = 4)+
  theme_bw()
  
ggsave("rawData.pdf", width = 20, height = 20, units = "cm")

data_nophage <- merged_data %>%
  filter(Treatment == "alone")  %>%
  filter(Sample_line == "Phage absent" | Sample_line == "Phage present") %>%
  filter(Time > 0) %>%
  filter(Time <= 1000)

set.seed(140520)

  
#y <-  (K*y0) / (y0 +(K-y0)  * exp(-mumax*time)) 

p     <- c(y0 = 0.1, mumax = 0.006, K = 0.6)
lower <- c(y0 = 0.08, mumax = 0.003, K = 0.5)  
upper <- c(y0 = 0.14, mumax = 0.007, K = 0.8)  
  
  
growth_model <- all_growthmodels(
  OD~grow_logistic(Time, parms) | SampleID + Sample_line, 
  data = data_nophage, 
  p = p, lower = lower, upper = upper, log = "y", ncores = 1) 

pp <- coef(growth_model)
  
growth_model2 <- all_growthmodels(
  OD~grow_logistic(Time, parms) | SampleID + Sample_line, 
  data = data_nophage, 
  p = pp, lower = lower, upper = upper, log = "y", ncores = 1)  
  
coef(growth_model2)  


par(mfrow = c(3,2))  
plot(growth_model2, log = "y")
  


p     <- c(y0 = 0.1, mumax = 0.006, K = 0.6, lambda =3)
lower <- c(y0 = 0.08, mumax = 0.003, K = 0.5, lambda = 1)  
upper <- c(y0 = 0.14, mumax = 0.007, K = 0.8, lambda = 10)  

  
growth_model_gompertz <- all_growthmodels(
  OD~grow_gompertz3(Time, parms) | SampleID + Sample_line, 
  data = data_nophage, 
  p = p, lower = lower, upper = upper, log = "y", ncores = 1)   

pp2 <- coef(growth_model_gompertz)

growth_model2_gompertz  <- all_growthmodels(
  OD~grow_gompertz2(Time, parms) | SampleID + Sample_line, 
  data = data_nophage, 
  p = pp2, lower = lower, upper = upper, log = "y", ncores = 1)  

par(mfrow = c(3,2))  
plot(growth_model2_gompertz, log = "y")


results_gompertz <- as.data.frame(results(growth_model2_gompertz))

results_gompertz_plot <-  results_gompertz %>%
  select(Sample_line, mumax, K) %>%
  gather("parameter", "value", -Sample_line) %>%
  group_by(Sample_line, parameter) %>%
  summarise(mean_value = mean(value), sd_value = sd(value))


ggplot(results_gompertz_plot, aes(x = Sample_line, y = mean_value, color = Sample_line))+
  geom_point(size=3)+
  geom_errorbar(aes(ymin = mean_value - sd_value, ymax = mean_value + sd_value), width = .2)+
  scale_color_manual(values = c("grey25", "darkred"))+
  facet_wrap(.~parameter, scales = "free_y")+
  labs( x = "Sample_line", y = "Parameter values", title = "Gompertz nonlinear growth model")+
  theme_bw()+
  theme(legend.position = "none")


t.test(K~Sample_line, results_gompertz)

################################################


data_withphage <- merged_data %>%
  filter(Treatment == "with phage") %>%
  filter(Sample_line == "Phage absent" | Sample_line == "Phage present") %>%
  filter(Time > 0) %>%
  filter(Time <= 1000)


auc <- function(x) {sum(cumsum(x)/2)}

data_split <- split(data_nophage["OD"], list(data_nophage$SampleID))
many_auc   <- sapply(data_split, auc)


data_split2 <- split(data_withphage["OD"], list(data_withphage$SampleID))
many_auc2   <- sapply(data_split2, auc)


control_data <- cbind.data.frame(Sample = names(many_auc), AUC = many_auc)
phage_data  <- cbind.data.frame(SampleID = names(many_auc2), Sample = gsub("\\+.*", "", names(many_auc2)),
AUC = many_auc2)

attributes <- data_withphage[, c("SampleID", "Sample_line", "Phage_line", "Local_line")]
unique_attributes<- distinct(attributes)

infectivity_data <- phage_data %>%
  full_join(control_data, by = "Sample") %>%
  full_join(unique_attributes, by = "SampleID") %>%
  mutate(infectivity = ((AUC.y - AUC.x)/AUC.y)*100)
  

ggplot(infectivity_data, aes(y = infectivity, x = Local_line, color = Sample_line, shape = Phage_line))+
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed")+
  geom_point(size=5, position = position_dodge(width = 0.5), alpha = 0.5)+
  labs(x = "Local phage line", y = "Phage infectivity", title ="Estimation of phage virulence/efectivity")+
  scale_color_manual("Sample line", values = c("grey25", "darkred"))+
  scale_shape_manual("Phage line", values = c(16, 17))+
  theme_bw()
  




