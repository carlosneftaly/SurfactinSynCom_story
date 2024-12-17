#'--- 
#'Title: Growth curve on TSB 
#'Author: Carlos N. Lozano Andrade
#'Aim: 
#'Date: Sept, 15th, 2020 
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
library(growthcurver)

###### importing data #################### 
data.GC <- read.csv('Syncom_growthCurves/syncom_growth_curvs.csv',
                    header=T,
                    sep=';', 
                    dec=",") 
############# Adjust the table into R-readable format#################### 

data.GC_clean <- data.GC %>% pivot_longer(
  cols = -Time,
  names_to = c('Strain', 'Treatment', 'Rep'),
  names_sep = '_',
  values_to = "Abs"
)  %>%
  group_by(Time, Strain, Treatment) %>% summarise(
    OD600_mean = mean(Abs), OD600_sd = mean(Abs)
  )


data.GC_clean$Strain <- as.factor(data.GC_clean$Strain) 
data.GC_clean$Treatment <- as.factor(data.GC_clean$Treatment) 

################# Plot ################################################### 



legend_labels <- c( "*Pedobacter*",
                    "*Rhodococcus*",
                    "*Stenotrophomonas*",
                    "*Chryseobacterium*",
                    "*sfp*",
                    "*srfAC*",
                    "*WT*"
              )


ggplot(data.GC_clean, aes (x=Time, y = OD600_mean, colour = Strain)) + 
  geom_point() + facet_grid(~Treatment) + geom_smooth()  + 
  labs ( y = expression('OD'[600]), x = 'Time (h)') + theme_bw() +
  theme(
    legend.title=element_blank(),axis.text=element_text(size=18),
    axis.title=element_text(size=18,face="bold"),
    legend.text = element_markdown(size = 18),
    strip.text.x = element_text(size=20, face="bold")) + 
  scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                               "#87C55F" , "#9EB9F3", "#9EB9F4" ), labels = legend_labels)


###################### 0.1 TSB #################################

data.01TSB <- data.GC_clean %>% filter (Treatment == '0.1TSB')

ggplot(data.01TSB, aes (x=Time, y = OD600_mean, colour = Strain)) + 
  geom_point()  + geom_smooth()  + 
  labs ( y = expression('OD'[600]), x = 'Time (h)') + theme_bw() +
  theme(
    legend.title=element_blank(),axis.text=element_text(size=18),
    axis.title=element_text(size=18,face="bold"),
    legend.text = element_markdown(size = 18),
    strip.text.x = element_text(size=20, face="bold")) + 
  scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                               "#87C55F" , "#9EB9F3", "#7EB9F9" ), labels = legend_labels)

############ Growth Curver ################################## 

data.GC <- data.GC %>% filter(Time <= 1600)

gc_out <- SummarizeGrowthByPlate(data.GC,
                                 bg_correct = "min", plot_file = T)
head(gc_out)



dataSum_clean <- gc_out %>% separate(col = sample, c('Strain', 'Media','Rep'),
                                     sep = '_')


########### Growth rate (r) 0.1TSB ####################### 

dat.r.01TSB <- dataSum_clean %>% filter(Media == '0.1TSB')
  
dat.r.01TSB$Strain <- as.factor(dat.r.01TSB$Strain)

dat.r.01TSB$Strain <- factor(dat.r.01TSB$Strain , 
                                     levels = c('D764', 'D763',  'D757' ,'D749', 'WT', 'sfp', 'srfAC'),
                                     labels = c("Chryseobacterium", "Stenotrophomonas", 
                                                "Rhodococcus" ,"Pedobacter", 'WT', 'sfp', 'srfAC' ))


############ ANOVA over Growth r
mod.r <- lm(r ~ Strain, data = dat.r.01TSB)
    shapiro.test(mod.r$residuals)
    leveneTest(mod.r)
anova(mod.r)

out <- LSD.test(mod.r,"Strain", group=T)


ggplot(dat.r.01TSB, aes (x=Strain, y = r, colour = Strain)) + 
      geom_jitter(size = 4, position=position_jitter(0.2)) +
      stat_summary(fun = mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, color="black") + 
      labs ( y = bquote('Growth rate'~(min^-1)), x = 'Strains') + theme_thesis() +
     scale_color_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                           "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                           "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                           "#968EB9", "#BDDCBC", "#DE9D8B",
                           "#C0C0C0")) +  theme(legend.position = "none") 
  
####### Carrying capacity ########## 

ggplot(dat.r.01TSB, aes (x=Strain, y = k, colour = Strain)) + 
  ylim(0, 0.8) + geom_jitter(size = 4, position=position_jitter(0.2)) +
  stat_summary(fun = mean, fun.min=mean, fun.max=mean, geom="crossbar", width=0.5, color="black") + 
  labs ( y = bquote('Carrying capacity (k)'), x = 'Strains') + theme_thesis() +
  scale_color_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                              "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                              "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                              "#968EB9", "#BDDCBC", "#DE9D8B",
                              "#C0C0C0")) +  theme(legend.position = "none") 


############## Only syncon members ##################

dat.r.01TSB_syncom <- dat.r.01TSB %>% filter(Strain %in% c('Chryseobacterium',
                                                           'Stenotrophomonas',
                                                           'Pedobacter',
                                                           'Rhodococcus'))

ggplot(dat.r.01TSB_syncom, aes (x=Strain, y = r, fill = Strain)) + 
geom_point(size = 5, shape = 21, position=position_jitter(width=0.1)) + 
  ylim(c(0.003, 0.013)) +
  stat_summary(fun = mean, fun.min=mean, fun.max=mean, geom="crossbar", 
               width=0.3, color="black") + 
  labs ( y = bquote('Growth rate'~(min^-1)), x = 'Strains') +
  scale_fill_manual(values=color_lakota, labels = legend_labels) + 
  theme_thesis() + theme(legend.position = "none")



ggsave('growth_rate_SynCom.pdf', 
       path = here::here("Syncom_growthCurves/figs/"),
       width=75*3, height=93.5*3, units = "mm", dpi = 600) 




ggplot(dat.r.01TSB_syncom, aes (x=Strain, y = auc_l, fill = Strain)) + 
  geom_point(size = 5, shape = 21, position=position_jitter(width=0.1)) + 
  ylim(c(200, 800)) +
  stat_summary(fun = mean, fun.min=mean, fun.max=mean, geom="crossbar", 
               width=0.3, color="black") + 
  labs ( y = bquote('Area under growth curve'), x = 'Strains')  +
  scale_fill_manual(values=color_lakota, labels = legend_labels) + 
  theme_thesis() + theme(legend.position = "none")


ggsave('AUC_SynCom.pdf', 
       
       path = here::here("Syncom_growthCurves/figs/"),
       width=75*3, height=93.5*3, units = "mm", dpi = 600) 

mod.r_syncom <- lm(r ~ Strain, data = dat.r.01TSB_syncom)
shapiro.test(mod.r_syncom$residuals)
leveneTest(mod.r_syncom)
anova(mod.r_syncom)

out <- HSD.test(mod.r_syncom,"Strain", group=T)

mod.r_syncomK <- lm(auc_l ~ Strain, data = dat.r.01TSB_syncom)
shapiro.test(mod.r_syncomK$residuals)
leveneTest(mod.r_syncomK)
anova(mod.r_syncomK)

out2 <- HSD.test(mod.r_syncomK,"Strain", group=T)
