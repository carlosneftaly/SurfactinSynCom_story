#'--- 
#'Title: E.16 Syncom + Bacillus 
#'Author: Carlos N. Lozano Andrade
#'Aim: 
#'Date: Oct, 27th, 200
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
library(MESS)

#### Importing ################### 

dat_E16_syncom <- read.csv('E.16 Syncom full/2020.10.17 Syncom composition full.csv',
                           header = T, sep = ';')


dat_E16_syncom$Strain <- factor(dat_E16_syncom$Strain, 
                        levels = c('764', '763', '749', '757'),
                        labels = c("Chryseobacterium", "Stenotrophomonas", "Pedobacter",
                                   "Rhodococcus"))


legend_labels <- c( "*Chryseobacterium*","*Stenotrophomonas*",
                    "*Pedobacter*" , "*Rhodococcus*"
)



############### Carryng capacity by treatments#############
## Sum all the counts ()

dat_k <- dat_E16_syncom %>% 
          group_by(Time, Treatment, Replicate) %>% 
          summarize(Kcapacity = sum(CFU),
                    KC_mean = mean(CFU, na.rm = T),
                    KC_sd = sd(CFU, na.rm = T), 
                    KC_se = KC_sd/sqrt(KC_mean))
### plot ###################### 
  
  ggplot(dat_k, aes(x = Time, y = (KC_mean))) +
                geom_point() + facet_grid(~Treatment)


## Sum all the counts by Time, Strain ()

dat_k2 <- dat_E16_syncom %>% 
  group_by(Time, Treatment, Strain) %>% 
  summarize(Kcapacity = log10(sum(CFU)),
            KC_mean = log10(mean(CFU, na.rm = T)),
            KC_sd = log10(sd(CFU, na.rm = T)), 
            KC_se = KC_sd/sqrt(KC_mean))
### plot ###################### 

ggplot(dat_k2, aes(x = Time, y = (KC_mean), ymin = (KC_mean) - (KC_se),
                   ymax = (KC_mean) + (KC_se), colour = Strain)) +
  geom_line() +
  geom_ribbon()  + facet_grid(~Treatment) + ylim(0,10)




##################### Syncom composition######################### 

              ######## Control###########
dat.Ctrl <-  dat_E16_syncom %>% filter(Treatment == 'CONTROL') %>%
            group_by(Time, Strain) %>%
            summarise(MeanCo = mean(CFU, na.rm = T)) %>%
            mutate(Fraction = MeanCo/sum(MeanCo))
            
            Sctrl <- ggplot(dat.Ctrl, aes(x=Time, y=Fraction, fill=Strain)) + 
              geom_area(alpha=0.6 , size=.1, colour="white") +  
              labs(x=NULL, y='Proportion', title = "SynCom Baseline", size = 18) + 
              theme(legend.position = "none",
                    legend.title=element_blank(),
                    axis.text=element_text(size=18),
                    axis.title=element_text(size=18,face="bold"),
                    legend.text = element_markdown(size = 18),
                    strip.text.x = element_text(size=18, face="bold"),
                     
                    
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) + 
              scale_x_continuous(breaks = c(1,3, 6, 9, 12, 14), expand = c(0.01, 0)) +  
              scale_y_continuous(expand = c(0, 0.01)) +  
              guides(fill = guide_legend(reverse = TRUE)) + 
              
              scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                         "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                         "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                         "#968EB9", "#BDDCBC", "#DE9D8B",
                                         "#C0C0C0"), labels = legend_labels) + 
              theme_thesis()
            
            
            
            
            
            
            
            ggsave('Syncom_control.pdf',  width=6, height=6, units = "in", dpi = 600) 
            
            ###################WT##########################################
            dat.WT <-  dat_E16_syncom %>% filter(Treatment == 'WT') %>%
              group_by(Time, Strain) %>%
              summarise(MeanCo = mean(CFU, na.rm = T)) %>%
              mutate(Fraction = MeanCo/sum(MeanCo))
            
            Swt <- ggplot(dat.WT, aes(x=Time, y=Fraction, fill=Strain)) + 
              geom_area(alpha=0.6 , size=.1, colour="white") +  
              labs(x='Days', y='Proportion', title = "SynCom + WT",  size = 18) + 
              theme(legend.position = "none",
                    legend.title=element_blank(),
                    axis.text=element_text(size=18),
                    axis.title=element_text(size=18,face="bold"),
                    legend.text = element_markdown(size = 18),
                    strip.text.x = element_text(size=18, face="bold"),
                    
                    
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) + 
              scale_x_continuous(breaks = c(1,3, 6, 9, 12, 14), expand = c(0.01, 0)) +  
              scale_y_continuous(expand = c(0, 0.01)) +  
              guides(fill = guide_legend(reverse = TRUE)) + 
              
              scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                         "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                         "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                         "#968EB9", "#BDDCBC", "#DE9D8B",
                                         "#C0C0C0"), labels = legend_labels) + 
              theme_thesis()
            ggsave('syncom+WT.pdf',  width=6, height=6, units = "in", dpi = 600) 
            
            ###################SFP##########################################
            dat.SFP <-  dat_E16_syncom %>% filter(Treatment == 'SFP') %>%
              group_by(Time, Strain) %>%
              summarise(MeanCo = mean(CFU, na.rm = T)) %>%
              mutate(Fraction = MeanCo/sum(MeanCo))
            
            SFP <- ggplot(dat.SFP, aes(x=Time, y=Fraction, fill=Strain)) + 
              geom_area(alpha=0.6 , size=.1, colour="white") +  
              labs(x='Days', y=NULL, title = "SynCom + sfp",  size = 18) + 
              theme(legend.position = "none",
                    legend.title=element_blank(),
                    axis.text=element_text(size=18),
                    axis.title=element_text(size=18,face="bold"),
                    legend.text = element_markdown(size = 18),
                    strip.text.x = element_text(size=18, face="bold"),
                    
                    
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) + 
              scale_x_continuous(breaks = c(1,3, 6, 9, 12, 14), expand = c(0.01, 0)) +  
              scale_y_continuous(expand = c(0, 0.01)) +  
              guides(fill = guide_legend(reverse = TRUE)) + 
              
              scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                         "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                         "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                         "#968EB9", "#BDDCBC", "#DE9D8B",
                                         "#C0C0C0"), labels = legend_labels) + 
              theme_thesis()
            
            ggsave('syncom+SFP.pdf',  width=6, height=6, units = "in", dpi = 600) 
            
            
            ###################SRFAC##########################################
            dat.srfAC <-  dat_E16_syncom %>% filter(Treatment == 'srfAC') %>%
              group_by(Time, Strain) %>%
              summarise(MeanCo = mean(CFU, na.rm = T)) %>%
              mutate(Fraction = MeanCo/sum(MeanCo))
            
            srfAC <- ggplot(dat.srfAC, aes(x=Time, y=Fraction, fill=Strain)) + 
              geom_area(alpha=0.6 , size=.1, colour="white") +  
              labs(x='Days', y=NULL, title = "SynCom + srfAC",  size = 18) + 
              theme(legend.position = "none",
                    legend.title=element_blank(),
                    axis.text=element_text(size=18),
                    axis.title=element_text(size=18,face="bold"),
                    legend.text = element_markdown(size = 18),
                    strip.text.x = element_text(size=18, face="bold"),
                    
                    
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) + 
              scale_x_continuous(breaks = c(1,3, 6, 9, 12, 14), expand = c(0.01, 0)) +  
              scale_y_continuous(expand = c(0, 0.01)) +  
              guides(fill = guide_legend(reverse = TRUE)) + 
              
              scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                         "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                         "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                         "#968EB9", "#BDDCBC", "#DE9D8B",
                                         "#C0C0C0"), labels = legend_labels) + 
              theme_thesis()
            
            ggsave('syncom+srfac.pdf',  width=6, height=6, units = "in", dpi = 600) 
            
            
            ###################ppsC##########################################
            dat.ppsC <-  dat_E16_syncom %>% filter(Treatment == 'ppsC') %>%
              group_by(Time, Strain) %>%
              summarise(MeanCo = mean(CFU, na.rm = T)) %>%
              mutate(Fraction = MeanCo/sum(MeanCo))
            
            SppsC <- ggplot(dat.ppsC, aes(x=Time, y=Fraction, fill=Strain)) + 
              geom_area(alpha=0.6 , size=.1, colour="white") +  
              labs(x='Days', y=NULL, title = "SynCom + ppsC", size = 18)+ 
              theme(legend.position = "none",
                    legend.title=element_blank(),
                    axis.text=element_text(size=18),
                    axis.title=element_text(size=18,face="bold"),
                    legend.text = element_markdown(size = 18),
                    strip.text.x = element_text(size=18, face="bold"),
                    
                    
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) + 
              scale_x_continuous(breaks = c(1,3, 6, 9, 12, 14), expand = c(0.01, 0)) +  
              scale_y_continuous(expand = c(0, 0.01)) +  
              guides(fill = guide_legend(reverse = TRUE)) + 
              
              scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                         "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                         "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                         "#968EB9", "#BDDCBC", "#DE9D8B",
                                         "#C0C0C0"), labels = legend_labels) + 
              theme_thesis()
            
          
            ggsave('syncom+ppsc.pdf',  width=6, height=6, units = "in", dpi = 600) 
            
            
            #bACILLUS 
            dat_bac <- read.csv('E.16 Syncom full/2020.10.17 Syncom composition full - BACILLUS.csv', 
                                header = T, sep =';', dec = ',')
            
            dat_bac$Treatment <- as.factor(dat_bac$Treatment)
            dat_bac$Strain <- as.factor(dat_bac$Strain)
            
            dat_bac$Treatment <- factor(dat_bac$Treatment, levels = c("WT", 
                                                                      "SFP",
                                                                      "SRFAC", 
                                                                      "PPSC"))
            
            ## Lines WT and sfp 
           # bac <- c('WT', 'sfp')
            legend_line <- c( "*WT*","*sfp*", "*srfAC*", "*ppsC*" )
            dat_bacS <- dat_bac %>% 
              
              group_by(Day, Treatment) %>%
              summarise(MeanCo = mean(log10(Count), na.rm = T), 
                        sD = sd(log10(Count), na.rm = T))
            
            bacFC <- ggplot(dat_bacS, aes(Day, MeanCo, fill = Treatment)) + 
              geom_line( size=1) + ylim(c(2,8)) +
              geom_ribbon(aes(ymin = MeanCo - sD, ymax = MeanCo+sD), alpha = 0.8) +
              labs (y = 'Log (CFU/g)', x = 'Time (days)', size = 18)   + 
              theme(legend.position = "right",
                    legend.title=element_blank(),
                    axis.text=element_text(size=18),
                    axis.title=element_text(size=18,face="bold"),
                    legend.text = element_markdown(size = 18),
                    strip.text.x = element_text(size=18, face="bold"),
                    
                    
                    panel.border = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()) + 
              scale_x_continuous(breaks = c(1,3, 6, 9, 12, 14), expand = c(0.01, 0)) +  
              guides(fill = guide_legend(reverse = TRUE)) + 
              
              scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                         "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                         "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                         "#968EB9", "#BDDCBC", "#DE9D8B",
                                         "#C0C0C0"), labels = legend_line) + 
              theme_thesis()
            
            
            
            ggsave('syncom_Bacillus.pdf',  width=8, height=6, units = "in", dpi = 300) 
            
            ####### Bacillus alone ############## 
            
            bac_alone <- read.csv('E.16 Syncom full/bacillus_alone.csv', header = T, 
                                  sep =';', dec=',')
            
            
            bac_alone$Treatment <- factor(bac_alone$Treatment, levels = c("WT", 
                                                                       "SFP",
                                                                       "SRFAC", 
                                                                       "PPS"))
            
            
            legend_line <- c( "*WT*","*sfp*", "*srfAC*", "*ppsC*" )
            bac_aloneS <- bac_alone %>% 
              
              group_by(Day, Treatment) %>%
              summarise(MeanCo = mean(log10(Count.1), na.rm = T), 
                        sD = sd(log10(Count.1), na.rm = T))
            
            bac_aloneplot <- ggplot(bac_aloneS, aes(Day, MeanCo, fill = Treatment)) + 
              geom_line( size=1) + 
              geom_ribbon(aes(ymin = MeanCo - sD, ymax = MeanCo+sD), alpha = 0.8) +
              labs (y = 'Log (CFU/g)', x = 'Time (days)', size = 18)   +
              theme(legend.title=element_blank(),
                    axis.text=element_text(size=18),
                    axis.title=element_text(size=18,face="bold"),
                    legend.text = element_markdown(size =15),
                    strip.text.x = element_text(size=18, face="bold"),axis.ticks.y = element_blank()) +
              ylim(1,8) + 
              scale_x_continuous(breaks = c(1,3, 6, 9, 12)) + 
              
              scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                         "#87C55F"), labels = legend_line) + theme_thesis()
            
            ggsave('syncom_Bacillus_alone.png',  width=8, height=6, units = "in", dpi = 300) 
            
            
            ## Arraging 
            
            
            
            
            layout <- c(
              
              area(1,3,1,3),
              area(1, 2, 1, 2)
            )
            
            # Show the layout to make sure it looks as it should
            plot(layout)
            
            
            areas <- c(area(1, 2, 1,5),
                       area(2, 1), area (2,2), 
                       area(2, 3), area(2,4), area(2,5))
            
            
            
            # Show the layout to make sure it looks as it should
            plot(areas)
            
            
            areas <- c(area(1, 1, 1,4))
            
            
            
            # Show the layout to make sure it looks as it should
            plot(areas)
            
            
            
            
            
           bagecoSynCom <- Sctrl + Swt +  SFP + srfAC + plot_layout(ncol = 4 ) 
            
           ggsave('syncom_bageco.pdf',  width=80, height=25, units = "cm", dpi = 300) 
            bacFC + Sctrl + Swt +  SFP + srfAC + SppsC+ plot_layout(design = areas)
            
            
            (  bacFC +  plot_layout(ncol = 1 ) ) / (Sctrl + Swt +  Ssfp)
            
            
            syncom_Bacillus <- ( plot_spacer()  +  
                                   bacFC  ) / (Sctrl + Swt +  Ssfp)
            
            
            
            
            ggsave('2021.06.13 syncom_Bacillus_full.png',  width=16, height=8, units = "in", dpi = 300) 
            
            
            areas2 <- c(area(1, 1, 1,5))
            plot(areas2)
            
            Sctrl + Swt +  SFP + srfAC + SppsC + plot_layout(design = areas2)
            
            
############################### AUC ############################################################
            
          ### NEsted dataframe
           
            
            syncom_growth <-dat_E16_syncom %>% 
              group_by(Treatment, Strain,Replicate) %>% nest() %>%
              mutate(
                growth_auc = map_dbl( 
                        data, ~MESS::auc(.x$Time, .x$CFU))
              )

            ### Summarizing the growth curve data 
            
            syncom_growth_sum <- syncom_growth %>% 
                        group_by(Treatment, Strain) %>%
                          summarize( 
                            syn_mean = mean(growth_auc),
                            syn_median = median(growth_auc),
                            syn_sd = sd(growth_auc), 
                            syn_err = (syn_sd/syn_mean)*100
                          )
            
            ggplot(syncom_growth_sum, aes(x = Treatment, y = log10(syn_mean ), colour = Strain )) + 
              geom_point() + geom_pointrange(aes(ymin = log10(syn_mean ) - log(syn_sd), ymax = log10(syn_mean ) + log(syn_sd) ))
            