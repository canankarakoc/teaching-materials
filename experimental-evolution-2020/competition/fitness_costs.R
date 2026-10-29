# Practical training UL 
# 12.05.20
# Karakoc
# Competition experiments/selection coefficient/fitness

# Load libraries
library(tidyverse) # we use dplyr and ggplot2

# Set working directory
# run from the folder this script is in

# Read data
comp_data <- read.table("competition_data.csv", header = T, sep = ";", dec = ",")

head(comp_data)
str(comp_data)

# Data manipulation
# For more https://dplyr.tidyverse.org/

new_comp_data <- comp_data %>% 
  filter(Dilution_fold == 10000 | Dilution_fold == 1000000) %>% # remove the highest dilutions
  mutate(Red_col_ml = Red_colonies * Dilution_fold * 100)  %>% # calculate CFU per mL
  mutate(Pink_col_ml = Pink_colonies * Dilution_fold * 100) %>% # calculate CFU per mL
  group_by(Time, SampleID, Sample_line) %>% # group the data to calculate the mean of replicate plates
  summarize(ancestral= mean(Red_col_ml), # calculate the mean of replicate plates
            evolved  = mean(Pink_col_ml)) %>% # calculate the mean of replicate plates
  mutate(total_biomass = ancestral + evolved) %>% # calculate total biomass
  gather("type", "cfu", -Time, -SampleID, -Sample_line) # create a longitudinal table for ggplot


# Plot 
ggplot(new_comp_data, aes(x = type, y = cfu, color = type))+
  geom_boxplot()+
  geom_point()+
  labs(x = "Colony type", y = "CFU/mL", title = "Raw data")+
  scale_color_manual("Colony type", values = c("darkred", "pink", "blue"))+
  facet_grid(Time~Sample_line, scales = "free")+
  theme_bw()
  
ggsave("rawData.pdf", width = 20, height = 15, units = "cm")


# Fitness costs
# Calculating selection coefficient 

fitness_data <- new_comp_data %>% 
  spread(type, cfu) %>% # create wide data table from long data table
  select(-total_biomass) %>% # remove the column total_biomass
  pivot_wider(names_from = "Time", 
              values_from = c("ancestral", "evolved")) %>% # make the data table even wider for calculations using Time variable
  mutate(selection_coef = (log(evolved_t1/ancestral_t1) - # calculate selection coefficient
                           log(evolved_t0/ancestral_t0))) %>%
  mutate(fitness = (selection_coef+1)) #calculate fitness
  
# Plot 
ggplot(fitness_data, aes(x = Sample_line, y = selection_coef, color = Sample_line))+
  geom_hline(yintercept=0, color="grey", linetype="dashed")+
  geom_boxplot()+
  geom_point()+
  labs(x = "Sample_line", y = "Selection coefficient", title = "s>0 means that the evolved have higher fitness")+
  scale_color_manual("Sample line", values = c("grey25", "darkred"))+
  theme_bw()

ggsave("rselection_coef.pdf", width = 15, height = 10, units = "cm")

# Plot 
ggplot(fitness_data, aes(x = Sample_line, y = fitness, color = Sample_line))+
  geom_hline(yintercept=0, color="grey", linetype="dashed")+
  geom_boxplot()+
  geom_point()+
  labs(x = "Sample line", y = "Fitness", title = "w>1 means increased fitness")+
  scale_color_manual("Sample line", values = c("grey25", "darkred"))+
  theme_bw()

ggsave("fitness.pdf", width = 15, height = 10, units = "cm")


# Alternative visualization with mean and SD

alternative_data <- fitness_data  %>% 
  select(SampleID, Sample_line, selection_coef, fitness) %>% 
  group_by(Sample_line) %>% # the grouping variable
  summarise(mean_selection_coef = mean(selection_coef), # calculates the mean selection coeffcient 
            sd_selection_coef = sd(selection_coef), # calculates the standard deviation of each group
            mean_fitness = mean(fitness), # calculates the mean fitness
            sd_fitness = sd(selection_coef)) # calculates the standard deviation of each group


# Plot 
ggplot(alternative_data, aes(x = Sample_line, y = mean_selection_coef, color = Sample_line))+
  geom_hline(yintercept=0, color="grey", linetype="dashed")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin = mean_selection_coef - sd_selection_coef, ymax = mean_selection_coef + sd_selection_coef), width=0.2)+
  labs(x = "Sample line", y = "Fitness", title = "w>1 means increased fitness")+
  scale_color_manual("Sample line", values = c("grey25", "darkred"))+
  theme_bw()+
  theme(legend.position = "none") #removes the legend

ggsave("fitness_alternative.pdf", width = 15, height = 10, units = "cm")

  
# Plot 
ggplot(alternative_data, aes(x = Sample_line, y = mean_fitness, color = Sample_line))+
  geom_hline(yintercept=0, color="grey", linetype="dashed")+
  geom_point(size=3)+
  geom_errorbar(aes(ymin = mean_fitness - sd_fitness, ymax = mean_fitness + sd_fitness), width=0.2)+
  labs(x = "Sample line", y = "Fitness", title = "w>1 means increased fitness")+
  scale_color_manual("Sample line", values = c("grey25", "darkred"))+
  theme_bw()+
  theme(legend.position = "none") #removes the legend

ggsave("fitness_alternative.pdf", width = 15, height = 10, units = "cm")


# Are the sample means significantly different?
t.test(fitness~Sample_line, data=fitness_data)

