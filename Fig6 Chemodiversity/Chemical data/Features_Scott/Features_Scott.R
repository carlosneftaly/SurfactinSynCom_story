#'--- 
#'Title: Fig_Chemistry_Features_Count - Scott
#'Author: Carlos N. Lozano Andrade
#'Aim: 
#'Date: July 27 
#'
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
# Importing data################ 


data_features <- read.csv('Chemical data/Features_Scott/Feature_Count_Scott.csv', 
                          header = T, sep = ';')


ggplot(data_features, aes(x = Time, y = Count,color = as.factor(Treatment))) + 
  geom_line() + 
  theme(#legend.position = "none",
        legend.title=element_blank(),
        axis.text=element_text(size=18),
        axis.title=element_text(size=18,face="bold"),
        legend.text = element_markdown(size = 18),
        strip.text.x = element_text(size=18, face="bold"),
        
        
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +  
  
  guides(fill = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                             "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                             "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                             "#968EB9", "#BDDCBC", "#DE9D8B",
                             "#C0C0C0")) + theme_thesis() + 
  scale_x_continuous(breaks = c(1,3,6,9))
ggsave('Features_Scott.pdf',  width=9, height=6, units = "in", dpi = 600) 

  