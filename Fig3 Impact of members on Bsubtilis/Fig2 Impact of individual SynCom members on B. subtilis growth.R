#'--- 
#'Title: E.21. Invtiro: Growth curves  
#'Author: Carlos N. Lozano Andrade
#'Aim: 
#'Date: Dec, 13th, 2021
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

#### Importing ################### 

dat.growth_gfp <- read.csv('Growth_curve/E20. growthcurve_isolates_D749_gfp.csv',
                           header = T, sep = ';', dec = ',')




plot(dat.growth_gfp$Time, dat.growth_gfp$srfAC_D1_R1)


trapz(dat.growth_gfp$Time, dat.growth_gfp$srfAC_D1_R1)




dat.growth_gfp <- dat.growth_gfp %>% filter(Time < 1800)


summG <- function(x) {trapz(dat.growth_gfp$Time,x)}
auc_trap <- sapply(dat.growth_gfp[2:ncol(dat.growth_gfp)], summG)

par(mfcol = c(8,11))
par(mar = c(0.25,0.25,0.25,0.25))
summGplot <- function(x) {
  
  
  
  plot(dat.growth_gfp$Time,x)
  
  
  
}
sapply(dat.growth_gfp[2:ncol(dat.growth_gfp)], summGplot)



df_auc_trap <- enframe(auc_trap)

df_auc_trap_clean <- df_auc_trap %>% separate(col = name, c('Strain', 'Dilution', 
                                                            'Rep'), sep = '_')


dat.sum.aud <- df_auc_trap_clean  %>% 
  group_by(Strain, Dilution) %>% summarise(AUC_mean = mean (value))




df_auc_trap_clean_mono <- df_auc_trap_clean %>% filter(Dilution == 'MONO')

dataMono_sum <- df_auc_trap_clean_mono %>% group_by(Strain) %>% summarise(AUC_mean = mean (value))


############## Relative to mono-culture at 10 , ######################### 

data_AUC_relative <- df_auc_trap_clean %>% 
  group_by(Strain, Dilution, Rep) %>%
  mutate(auc_fold = 
           case_when(
             Strain == 'WT' ~ value/22275327,
             Strain == 'sfp' ~ value/73560013,
             Strain == 'srfAC' ~ value/35575268
           )
  )  




data_AUC_relative$Strain <- factor(data_AUC_relative$Strain, 
                                   levels = c('WT', 'sfp',  'srfAC'))


data_AUC_relative %>% filter(Dilution == 'D1')

ggplot(data_AUC_relative %>% filter(Dilution == 'MONO'), aes(x= Strain, y = auc_fold, fill = Strain)) + 
  geom_boxplot() + geom_jitter()  + 
  labs(y = 'AUDC', title = 'MONO')  +
  theme(legend.position = "none",
        legend.title=element_blank(),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"),
        legend.text = element_markdown(size = 22),
        strip.text.x = element_text(size=22, face="bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +  
  guides(fill = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                             "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                             "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                             "#968EB9", "#BDDCBC", "#DE9D8B",
                             "#C0C0C0")) 



ggplot(data_AUC_relative %>% filter(auc_fold < 3) , aes(x= Strain, y = auc_fold, fill = Strain)) + 
  geom_boxplot() + geom_jitter()  + facet_grid(~Dilution) + 
  labs(y = 'AUDC_fold')  +
  theme(legend.position = "none",
        legend.title=element_blank(),
        axis.text=element_text(size=16),
        axis.title=element_text(size=22,face="bold"),
        legend.text = element_markdown(size = 22),
        strip.text.x = element_text(size=22, face="bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +  
  guides(fill = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                             "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                             "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                             "#968EB9", "#BDDCBC", "#DE9D8B",
                             "#C0C0C0")) 




dat.spent_OD600_s10 <- dat.spent_OD600 %>% select(contains(c('Time', 'S10', 'MONO')))


p <- lapply(
  colnames(dat.spent_OD600_s10)[2:ncol(dat.spent_OD600_s10)],
  function(col) ggplot(dat.spent_OD600_s10, aes_string(x = 'Time', y = col)) +
    geom_point()  +
    theme_bw()) 


require(cowplot)
margin <- theme(
  plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), 'cm'),
  plot.title = element_blank(),
  plot.subtitle = element_blank(),
  plot.caption = element_blank())


plots <- lapply(p, '+', margin)


cairo_pdf('full_plot.pdf', width =24, height = 8)
plot_grid(
  do.call(plot_grid, c(plots, ncol = 12, nrow = 4)),
  ncol = 1)

dev.off()




############## Relative to mono-culture at 10 , ######################### 

data_AUC_S10_relative <- dataAUC_s10 %>% 
  group_by(Strain, SynMember, Rep) %>%
  mutate(auc_fold = 
           case_when(
             Strain == 'WT' ~ value/385.8283,
             Strain == 'SFP' ~ value/1000.3075,
             Strain == 'srfAC' ~ value/441.0192
           )
  )  




data_AUC_S10_relative$SynMember <- factor(data_AUC_S10_relative$SynMember, 
                                          levels = c('D764', 'D763',  'D757' ,'D749'),
                                          labels = c("Chryseobacterium", "Stenotrophomonas", 
                                                     "Rhodococcus" ,"Pedobacter"))

data_AUC_S10_relative$Strain <- factor(data_AUC_S10_relative$Strain, 
                                       levels = c('WT', 'SFP',  'srfAC'))




ggplot(data_AUC_S10_relative, aes(x= Strain, y = auc_fold, fill = Strain)) + 
  geom_boxplot() + geom_jitter() + facet_grid(~SynMember)  + ylim(0.2, 1.5) + 
  labs(y = 'AUDC relative to Monoculture', title = 'OD600')  +
  theme(legend.position = "none",
        legend.title=element_blank(),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"),
        legend.text = element_markdown(size = 22),
        strip.text.x = element_text(size=22, face="bold"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())  +  
  guides(fill = guide_legend(reverse = TRUE)) + 
  scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                             "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                             "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                             "#968EB9", "#BDDCBC", "#DE9D8B",
                             "#C0C0C0")) + geom_hline(yintercept=1, size=1.5) 




ggsave('AUDC_OD600_trapezoid.png',  width=12, height=7, units = "in", dpi = 300) 






### Heat maps ######################## 

ggplot(data_AUC_S10_relative, aes(x=Strain, y=SynMember)) + 
  geom_tile(colour="black", size=0.25, aes(fill=auc_fold)) +
  scale_fill_gradient(low="#ecefb7", high="#3a81b5", na.value = "white") +
  theme(axis.text.y = element_text(size = 8))




## plotting all the columns same time 

summG <- function(x) {SummarizeGrowth(dat.spent$time,x)}
lapply(dat.spent[2:ncol(dat.spent)], summG)


models.all <- lapply(dat.spent[2:ncol(dat.spent)], function(x) SummarizeGrowth(dat.spent$time, x))

df.predicted.plate <- data.frame(time = dat.spent$time)
for (i in names(dat.spent[2:ncol(dat.spent)])) 
{df.predicted.plate[[i]] <- predict(models.all[[i]]$model)}


melt1 <- melt(dat.spent, id.vars = "time", variable.name = "sample", value.name = "GFP")
melt2 <- melt(df.predicted.plate, id.vars = "time", variable.name = "sample", value.name = "pred.od")
df.final <- cbind(melt1, pred.od=melt2[,3])

ggplot(df.final, aes(x=time, y=GFP)) + 
  geom_point(aes(), alpha=0.5) + geom_line(aes(y=pred.od), color="red") + 
  facet_wrap(~sample, ncol = 12) + theme_bw()



######## Tile plot ############ 


### Reduction (%)###########3
# auc_red = (growtMono - growthColculture)/GrowthMono)*100
dataMono_sum <- df_auc_trap_clean_mono %>% group_by(Strain) %>% summarise(AUC_mean = median (value))


data_AUC_red<- df_auc_trap_clean %>% 
  group_by(Strain, Dilution, Rep) %>%
  mutate(auc_red = 
           case_when(
             Strain == 'WT' ~ ((21638785 - value)/(21638785)),
             Strain == 'sfp' ~ ((76547900 - value)/(76547900)),
             Strain == 'srfAC' ~ ((33730250 - value)/(33730250))
           )
  ) 

## eWith medians
data_tile <- data_AUC_red %>% 
  group_by(Strain, Dilution) %>% 
  summarise(AUC_red = median (auc_red))


data_tile



ggplot(data_tile, aes(x=Strain, y=SynCom)) + 
  geom_tile(colour="grey50", size=0.25, aes(fill=AUC_red)) +
  scale_fill_gradient(low="#3388bd", high="#e05753", na.value = "white",
                      breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), 
                      labels = c("0.0", "0.25", "0.5", "0.75", "1.0")) +
  theme(axis.text.y = element_text(size = 12)) + theme_minimal()





dat.sum.aud <- df_auc_trap_clean  %>% 
  group_by(Strain, Dilution) %>% summarise(AUC_mean = mean (value))
