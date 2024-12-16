#'--- 
#'Title: Fig_Syncom Composition
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
library(ComplexHeatmap) #install_github("jokergoo/ComplexHeatmap")
library(pheatmap)
library("FactoMineR")
library("factoextra")
library(vegan)
library(glue)
# Importing data################ 

data_compostion <- read.csv('Composition_PCA/Syncom_composition_PCA.csv',
                            header = T, sep = ';')

comp.pca <- PCA(data_compostion[,4:7], graph = T)


fviz_pca_ind(comp.pca,
             geom.ind = "point", # show points only (nbut not "text")
             #fill.ind = c(as.factor(data_compostion$Treatment)),
             col.ind = c(as.factor(data_compostion$Treatment)), 
             #habillage = as.factor(data_compostion$Treatment),# color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07", 'black', 'green', 'yellow'),
             addEllipses = F, # Concentration ellipses
             legend.title = "Groups"
             
)



comp.pca2 <- prcomp(data_compostion[,4:7],  scale = TRUE)
basic_plot <- fviz_pca_ind(comp.pca2, label="none")

ggplot(cbind(basic_plot$data,data_compostion[,c("Time","Treatment")]),
       aes(x=x,y=y,col=Treatment,shape=as.factor(Time))) + geom_point() + theme_bw() + xlim(-2,6)






############## NMDS ############ 
set.seed(19880604)
nmds <- metaMDS(as.matrix(data_compostion[,4:7]), distance = "bray")
nmds
plot(nmds)
data.scores = as.data.frame(scores(nmds)$sites)

#add columns to data frame 
data.scores$Time = as.factor(data_compostion$Time)
data.scores$Treatment = data_compostion$Treatment


ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Treatment, colour = Time)) + 
  
  scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                             "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                             "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                             "#968EB9", "#BDDCBC", "#DE9D8B",
                             "#C0C0C0"))  + xlim(-1, .6) + 
  
  theme(
    
    legend.key = element_rect(fill = "#fbf7f1", colour = "#fbf7f1"),
    legend.key.size = unit(1.21, "cm"),
    legend.title = element_markdown(size=22, face="bold"),
    legend.text = element_markdown(size = 22),
    legend.margin = margin(c(5, 5, 5, 5)),
    axis.text.x = element_text(size=19, margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.text.y = element_text(size=19, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.title.x=element_text(size=25,face="bold", margin = margin(t = 10, r = 0, b = 0, l = 0)),
    axis.title.y=element_text(size=25,face="bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
    strip.background = element_rect(colour = "black", fill = "#F0F0F0", size = 2),
    strip.text.x = element_text(size=25, face="bold"),
    panel.spacing = unit(1, "lines"),
    panel.background = element_rect(fill = "#fbf7f1"),
    #panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
    panel.border = element_rect(colour = "black", fill=NA, size = 2), 
    plot.background = element_rect(fill = "transparent", colour = NA))


ggsave('nmds_composition2.pdf', width = 9, 
       height = 6, units = 'in', dpi = 600)

set.seed(123)
comp.permanova1 <- adonis2(data_compostion[,4:7] ~ as.factor(data_compostion$Time)*data_compostion$Treatment,
                          permutations = 9999, method="bray")
comp.permanova1
      
############ Removing day 1 ################## 
data_comp.d <- data_compostion %>% filter(Time != 1)

nmds <- metaMDS(as.matrix(data_comp.d[,4:7]), distance = "bray")
nmds
plot(nmds)
data.scores = as.data.frame(scores(nmds)$sites)

#add columns to data frame 
data.scores$Time = as.factor(data_comp.d$Time)
data.scores$Treatment = data_comp.d$Treatment


plot(data.scores[,1:2], col=factor(data.scores$Treatment))

ordispider(nmds, groups = data.scores$Treatment, col = 1:4)
ordiellipse(nmds, groups = data.scores$Treatment, col = 1:5)

ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Treatment, colour = Time)) + 
  
  scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                               "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                               "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                               "#968EB9", "#BDDCBC", "#DE9D8B",
                               "#C0C0C0")) + theme_thesis() + xlim(-1, .6)
set.seed(123)
comp.permanova1 <- adonis2(vegdist(data_compostion[,4:7]) ~ as.factor(data_compostion$Time)*data_compostion$Treatment,
                           permutations = 999, method="bray")
comp.permanova1

#REmoving control #################### 

     data_composition2 <- data_compostion %>% filter(Treatment != 'CONTROL')
      nmds2 <- metaMDS(as.matrix(data_composition2[,4:7]), distance = "bray")
      nmds2
      plot(nmds2)
      data.scores2 = as.data.frame(scores(nmds2)$sites)
      
      #add columns to data frame 
      data.scores2$Time = as.factor(data_composition2$Time)
      data.scores2$Treatment = data_composition2$Treatment
      
      
      ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(size = 4, aes( colour = Time, shape = Treatment)) + 
        
        scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                     "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                     "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                     "#968EB9", "#BDDCBC", "#DE9D8B",
                                     "#C0C0C0")) + xlim(-1, .4) + theme_poster() 
      
      ggsave('nmds.png', width = 25.2, height = 15.75, units = 'cm', dpi = 600, bg='#EAE8DC')
      
      
      ######## ANOSIM and Permanova ############## 
      data_composition2[,4:7]
    
      comp.permanova <- adonis2(vegdist(data_composition2[,4:7]) ~ as.factor(data_composition2$Time)*data_composition2$Treatment,
                               permutations = 999, method="bray")
      comp.permanova
      
      ## Day 1 
      
      data_day1 <- data_composition2 %>% filter(Time == 1)
      data_day1.permanova <- adonis2(data_day1[,4:7] ~ data_day1$Treatment,
                                permutations = 999, method="bray")
      data_day1.permanova
      
      
      nmds.day1 <- metaMDS(as.matrix(data_day1[,4:7]), distance = "bray")
      nmds.day1
      plot(nmds.day1)
      data.scores2 = as.data.frame(scores(nmds.day1)$sites)
      
      #add columns to data frame 
      data.scores2$Time = as.factor(data_day1$Time)
      data.scores2$Treatment = data_day1$Treatment
      
      
      ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(size = 4, aes( colour = Treatment)) + 
        
        scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                     "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                     "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                     "#968EB9", "#BDDCBC", "#DE9D8B",
                                     "#C0C0C0")) + theme_thesis()
      
      
      # Dy 3 
      data_day3 <- data_composition2 %>% filter(Time == 3)
      data_day3.permanova <- adonis2(data_day3[,4:7] ~ data_day3$Treatment,
                                     permutations = 1e4, method="bray")
      data_day3.permanova
      
      
      nmds.day3 <- metaMDS(as.matrix(data_day3[,4:7]), distance = "bray")
      nmds.day3
      plot(nmds.day3)
      data.scores2 = as.data.frame(scores(nmds.day3)$sites)
      
      #add columns to data frame 
      data.scores2$Time = as.factor(data_day1$Time)
      data.scores2$Treatment = data_day3$Treatment
      
      
      ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(size = 4, aes( colour = Treatment)) + 
        
        scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                     "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                     "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                     "#968EB9", "#BDDCBC", "#DE9D8B",
                                     "#C0C0C0")) + theme_thesis()
      
      
      
      
      
      
      
      
      # Day6
      data_day6 <- data_composition2 %>% filter(Time == 6)
      data_day6.permanova <- adonis2(data_day6[,4:7] ~ data_day6$Treatment,
                                     permutations = 999, method="bray")
      data_day6.permanova
      nmds.day6 <- metaMDS(as.matrix(data_day6[,4:7]), distance = "bray")
      nmds.day6
      plot(nmds.day6)
      data.scores2 = as.data.frame(scores(nmds.day6)$sites)
      
      #add columns to data frame 
      data.scores2$Time = as.factor(data_day1$Time)
      data.scores2$Treatment = data_day6$Treatment
      
      
      ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(size = 4, aes( colour = Treatment)) + 
        
        scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                     "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                     "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                     "#968EB9", "#BDDCBC", "#DE9D8B",
                                     "#C0C0C0")) + theme_thesis()
      # Day9
      data_day9 <- data_composition2 %>% filter(Time == 9)
      data_day9.permanova <- adonis2(data_day9[,4:7] ~ data_day9$Treatment,
                                     permutations = 999, method="bray")
      data_day9.permanova
      
      nmds.day9 <- metaMDS(as.matrix(data_day9[,4:7]), distance = "bray")
      nmds.day9
      plot(nmds.day9)
      data.scores2 = as.data.frame(scores(nmds.day9)$sites)
      
      #add columns to data frame 
      data.scores2$Time = as.factor(data_day1$Time)
      data.scores2$Treatment = data_day9$Treatment
      
      
      ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(size = 4, aes( colour = Treatment)) + 
        
        scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                     "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                     "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                     "#968EB9", "#BDDCBC", "#DE9D8B",
                                     "#C0C0C0")) + theme_thesis()
      
      
      
      
      # Day12
      data_day12 <- data_composition2 %>% filter(Time == 12)
      data_day12.permanova <- adonis2(data_day12[,4:7] ~ data_day12$Treatment,
                                     permutations = 999, method="bray")
      data_day12.permanova
      
      nmds.day12 <- metaMDS(as.matrix(data_day12[,4:7]), distance = "bray")
      nmds.day12
      plot(nmds.day12)
      data.scores2 = as.data.frame(scores(nmds.day12)$sites)
      
      #add columns to data frame 
      data.scores2$Time = as.factor(data_day1$Time)
      data.scores2$Treatment = data_day1$Treatment
      
      
      ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + 
        geom_point(size = 4, aes( colour = Treatment)) + 
        
        scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                     "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                     "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                     "#968EB9", "#BDDCBC", "#DE9D8B",
                                     "#C0C0C0")) + theme_thesis()
      
      
      
      
      
      # Day14
      data_day14 <- data_composition2 %>% filter(Time == 14)
      data_day14.permanova <- adonis2(data_day12[,4:7] ~ data_day14$Treatment,
                                      permutations = 999, method="bray")
      data_day14.permanova
      
      
      
      
      
      
      
      
      
      
######### PCOA ################## 
      dist.pcoa <- vegdist(log(data_compostion[,4:7]+1), method="bray") 
      pcoa <- cmdscale(dist.pcoa, k = 2, eig =T , add = T) 
      
      positions <- pcoa$points
      colnames(positions) <- c('pcoa1', 'pcoa2')
      
      percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)
      
      pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
      
      labels <- c(glue("PCo Axis 1 ({pretty_pe[1]}%)"),
                  glue("PCo Axis 2 ({pretty_pe[2]}%)"))
      
      positions %>%
        as_tibble() %>%
        ggplot(aes(x=pcoa1, y=pcoa2)) +
        geom_point(size = 4, aes( colour = data_compostion$Treatment, 
                                  shape = as.factor(data_compostion$Time))) +
        labs(x=labels[1], y=labels[2])
      
      
      
      
     
      tibble(pe = cumsum(percent_explained),
             axis = 1:length(percent_explained)) %>%
        ggplot(aes(x=axis, y=pe)) +
        geom_line() +
        coord_cartesian(xlim = c(1, 10), ylim=c(0, 100)) +
        scale_x_continuous(breaks=1:10)
      
      
      
    
    
    
    library(ape)
    PCOA <- pcoa(dist.pcoa)
    PCOA$values
    biplot.pcoa(PCOA)
    