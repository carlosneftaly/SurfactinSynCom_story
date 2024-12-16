#'--- 
#'Title: Fig_Syncom induces Surfactin production
#'Author: Carlos N. Lozano Andrade
#'Aim: 
#'Date: June 6th
#'--


#### Libraries #################
# The next libraries MUST be installed (the first time) or loaded before start running the script
library(tidyverse)
library(scales)
library(agricolae)
library(car)
library(ggpval)
library(ggsci)
library(patchwork)
library(ggtext)
library(ComplexHeatmap)
library(pheatmap)
library("FactoMineR")
library("factoextra")
library(vegan)
library(MetBrewer)
# Importing data################ 
color_lakota <- met.brewer("Lakota", 8, override.order=TRUE)[c(8,7,5,4,2,1)]  


data_induction<- read.csv('FigX_Syncom induces Surfactin_/Surfactin quant_Carlos 051120.csv',
                            header = T, sep = ';', dec =',')

ggplot(data_induction, aes(x = Condition, y = Concentration))  + 
  geom_boxplot(aes(fill = Condition)) + 
  geom_point(size = 3, shape = 21, fill="grey20",
             position=position_jitter(width=0.2, height=0.1)) +

  guides(fill = guide_legend(reverse = TRUE)) + 
  ylim(30, 80)  +
  scale_fill_manual(values=color_lakota) + 
  theme_thesis() + theme (legend.position = "none")

ggsave('Induction_boxplot.pdf',  
       path = here::here("FigX_Syncom induces Surfactin_/figs/"),
       width=5, height=6, units = "in", dpi = 600) 



######## Anova ################## 

mod.induces <- lm(Concentration~Condition, data = data_induction)

shapiro.test(mod.induces$residuals)
leveneTest(Concentration~Condition, data = data_induction)
anova(mod.induces)

t.test(Concentration~Condition, data = data_induction, var.equal = T)

kruskal.test(Concentration~Condition, data = data_induction)

######### Supernantant Caja thesis ################## 

data_hplc_CFS <- read.csv('FigX_Syncom induces Surfactin_/SM_induction_HPLC_supernantant.csv',
                          header = T, sep = ';', dec = ',')




data_hplc_CFS$Treatment <- as.factor(data_hplc_CFS$Treatment)



data_sum <- data_hplc_CFS %>% group_by(Treatment) %>%
              summarise(Fold_mean = mean(Area_foldM),
                        Fold_median = median(Area_foldM))




ggplot() + 
  geom_bar(data = data_sum, aes(x = fct_rev(Treatment), y = Fold_median, fill = fct_rev(Treatment)),
           stat="identity", position = position_dodge())   +
  scale_fill_manual(values=color_lakota) + 
  theme_thesis() + theme(legend.position = "none") + 

  geom_point(data = data_hplc_CFS,aes(x = fct_rev(Treatment), y = Area_foldM, fill = fct_rev(Treatment)),
             size = 5, shape = 21, fill="grey20",
             position=position_jitter(width=0.2, height=0.1)) + coord_flip() +
  geom_hline(yintercept = 1,  color = "black", size=1.5) 

  
  ggsave('CFS_CajaTHesis.pdf', 
         path = here::here("FigX_Syncom induces Surfactin_/figs/"), 
         width=9, height=6, units = "in", dpi = 600) 
  
  