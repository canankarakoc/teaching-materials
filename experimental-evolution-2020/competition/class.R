# Praktical training UL 
# 12.05.20
# Karakoc
# Competition experiments/selection coefficients/fitness

# Load libraries 
library(tidyverse)  # we use dplyr and ggplot2 

# Set working directory 
# run from the folder this script is in

# Read data 
comp_data <- read.table("competition_data.csv", header = T, sep = ";", dec = ",")

head(comp_data)
str(comp_data)

# Data manipulation 

new_comp_data <- comp_data %>%
  filter(Dilution_fold == 10000 | Dilution_fold == 1000000) %>%
  mutate(Red_col_ml = Red_colonies * Dilution_fold * 100) %>%
  mutate(Pink_col_ml = Pink_colonies * Dilution_fold * 100)  %>%
  group_by(Time, SampleID, Sample_line)  %>%
  summarise(ancestral = mean(Red_col_ml), evolved = mean(Pink_col_ml))  %>%
  mutate(total_biomass = ancestral + evolved ) %>%
  gather("type","cfu", -Time, -SampleID, -Sample_line)

  
ggplot(new_comp_data, aes(x = type, y=cfu, color = type))+
  geom_boxplot()+
  geom_point()+
  labs(x="Colony type", y="CFU/mL", title = "Raw data")+
  scale_color_manual("Colony type", values = c("darkred", "pink", "blue"))+
  facet_grid(Time~Sample_line, scales="free")+
  theme_bw()
  
ggsave("rawData.pdf", width = 20, height = 15, units="cm")  
  
  
#Calculating selection coefficient 

fitness_data <- new_comp_data %>%
  spread (type, cfu) %>%
  select(-total_biomass) %>%
  pivot_wider(names_from = "Time", values_from = c("ancestral", "evolved")) %>%
  mutate(selection_coef = (log(evolved_t1/ancestral_t1)) - (log(evolved_t0/ancestral_t0)))%>%
  mutate(fitness = selection_coef + 1)
  
#Plot 
ggplot(fitness_data, aes(x = Sample_line, y = selection_coef, color = Sample_line))+
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed")+
  geom_boxplot()+
  geom_point()+
  labs(x = "Sample line", y = "Selection coefficient (s)", title ="s>0 means that evolved have hihger fitness")+
  scale_color_manual("Sample line", values = c("grey25", "darkred"))+
  theme_bw()


ggplot(fitness_data, aes(x = Sample_line, y = fitness, color = Sample_line))+
  geom_hline(yintercept = 1, color = "grey", linetype = "dashed")+
  geom_boxplot()+
  geom_point()+
  labs(x = "Sample line", y = "Fitness (w)", title ="w>1 means increased fitness")+
  scale_color_manual("Sample line", values = c("grey25", "darkred"))+
  theme_bw()

#Alternative method 

alternative_data <- fitness_data %>%
  select(SampleID, Sample_line, selection_coef , fitness) %>%
  group_by(Sample_line) %>%
  summarise(mean_selection_coef = mean(selection_coef), 
            sd_selection_coef = sd(selection_coef), 
            mean_fitness = mean(fitness), 
            sd_fitness = sd(fitness))
           


#Plot
ggplot(alternative_data, aes(x = Sample_line, y = mean_selection_coef, color = Sample_line))+
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed")+
  geom_errorbar(aes(ymin = mean_selection_coef - sd_selection_coef, ymax = mean_selection_coef + sd_selection_coef), width =0.2)+
  geom_point(size = 3)+
  labs(x = "Sample line", y = "Selection coefficient (s)", title ="s>0 means that evolved have hihger fitness")+
  scale_color_manual("Sample line", values = c("grey25", "darkred"))+
  theme_bw()+
  theme(legend.position = "none")

ggplot(alternative_data, aes(x = Sample_line, y = mean_fitness, color = Sample_line))+
geom_hline(yintercept = 1, color = "grey", linetype = "dashed")+
  geom_errorbar(aes(ymin = mean_fitness - sd_fitness, ymax = mean_fitness + sd_fitness ), width =0.2)+
  geom_point(size = 3)+
  labs(x = "Sample line", y = "Fitness (w)", title ="w>1 means increased fitness")+
  scale_color_manual("Sample line", values = c("grey25", "darkred"))+
  theme_bw()+
  theme(legend.position = "none")


t.test(fitness~Sample_line, data = fitness_data)



