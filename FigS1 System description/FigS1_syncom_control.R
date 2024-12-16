#'--- 
#'Title: Fig Sup1_ Syncom control 
#'Author: Carlos N. Lozano Andrade
#'Aim: Syncom without Bacillus inoculation /data previous 
#'Date: June, 1th, 2022
#'--

#### Libraries #################
# The next libraries MUST be installed (the first time) or loaded before start running the script
library(tidyverse)
library(scales)
library(ggpval)
library(ggsci)
library(patchwork)
library(ggtext)
library(growthcurver)
library(pracma)
library(svglite)
#### Importing ################### 

dat.syncom.C <- read.csv('Fig_Supp1_Syncom_Control/FigS1_syncom_control.csv',
                            header = T, sep = ';', dec = ',')
dat.Ctrl <-  dat.syncom.C %>% 
  group_by(Time, Strain) %>%
  summarise(MeanCo = mean(Count, na.rm = T)) %>%
  mutate(Fraction = MeanCo/sum(MeanCo))


dat.Ctrl$Strain <- factor(dat.Ctrl$Strain, 
                                levels = c('764', '763', '749', '757'),
                                labels = c("Chryseobacterium", "Stenotrophomonas", "Pedobacter",
                                           "Rhodococcus"))

legend_labels <- c( "*Chryseobacterium*","*S. indicatrix*",
                    "*Pedobacter*" , "*Rhodococcus*"
)

ggplot(dat.Ctrl, aes(x=Time, y=Fraction, fill=as.factor(Strain))) + 
  geom_area(alpha=0.6 , size=.1, colour="white") + 
  scale_x_continuous(breaks = c(1,3, 6, 9, 12, 14), expand = c(0.01, 0)) +  
  scale_y_continuous(expand = c(0, 0.01)) +
  guides(fill = guide_legend(reverse = TRUE)) + labs(x = 'Days') + 
  
  scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                             "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                             "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                             "#968EB9", "#BDDCBC", "#DE9D8B",
                             "#C0C0C0"), labels = legend_labels) +
  theme_thesis()


################# Final version ############################################# 

ggplot( dat.Ctrl, aes(x=Time, y=Fraction, fill=as.factor(Strain))) + 
  geom_area(alpha=0.7 , size=.2, colour="#C0C0C0") +  
  labs(x="Time (days)", y='Proportion') + 
  scale_x_continuous(breaks = c(1,3, 6, 9, 12, 14), expand = c(0.01, 0)) +  
  scale_y_continuous(expand = c(0, 0.01)) +  
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values=color_lakota, labels = legend_labels) + 
  theme_thesis() + theme(legend.position = "none")




ggsave('FigS1_Syncom_control_area.pdf', 
       path=here::here("Fig_Supp1_Syncom_Control" ,"figs"),
       width=6, height=8, units = "in", 
       dpi = 600, bg = "transparent") 

############# AS BArplot ################


dat.CtrlBar <-  dat.syncom.C %>% 
  group_by(Time, Strain) %>%
  summarise(MeanCo = mean(LogC, na.rm = T),
            StdCo = sd(LogC, na.rm = T))

dat.CtrlBar$Strain <- factor(dat.CtrlBar$Strain, 
                   levels = c('764', '763', '749', '757'),
                   labels = c("D764", "D763", "D749",
                              "D757"))

ggplot(data = dat.CtrlBar, aes(x = as.factor(Strain), y = MeanCo, fill = as.factor(Strain))) + 
  geom_bar(stat="identity", position = position_dodge()) + # adding bar plot
  geom_errorbar(aes(ymin=MeanCo-StdCo, ymax=MeanCo+StdCo),  
                width=.2, size = .8,                    # Width of the error bars
                position=position_dodge(.9)
  ) + facet_grid(~as.factor(Time)) + theme_thesis() + 
  
  scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                             "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                             "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                             "#968EB9", "#BDDCBC", "#DE9D8B",
                             "#C0C0C0")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  
  labs(x='Strain', y = 'Log10 (CFU/g)') 





########################### Final version #########################################

dat.CtrlBar <-  dat.syncom.C %>% 
  group_by(Time, Strain) %>%
  summarise(MeanCo = mean(LogC, na.rm = T),
            StdCo = sd(LogC, na.rm = T))

dat.CtrlBar$Strain <- factor(dat.CtrlBar$Strain, 
                             levels = c('764', '763', '749', '757'),
                             labels = c("D764", "D763", "D749",
                                        "D757"))

ggplot(data = dat.CtrlBar, aes(x = as.factor(Strain), y = MeanCo, fill = as.factor(Strain))) + 
  geom_bar(stat="identity", position = position_dodge()) + # adding bar plot
  geom_errorbar(aes(ymin=MeanCo-StdCo, ymax=MeanCo+StdCo),  
                width=.2, size = .8, colour = 'grey80',           # Width of the error bars
                position=position_dodge(.9)) + 
  facet_grid(~as.factor(Time))  +
  scale_fill_manual(values=color_lakota, labels = legend_labels) + 
  theme_thesis() + theme(legend.position = "none",
                         axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) + 
  labs(x='Strain', y = 'Log10 (CFU/g)') + 
 
  scale_y_continuous(expand = c(0, 0.01)) +  
  guides(fill = guide_legend(reverse = TRUE))

ggsave('FigS1_Syncom_control.pdf', 
       path=here::here("Fig_Supp1_Syncom_Control" ,"figs"),
       width=12, height=8, units = "in", dpi = 600) 




  
