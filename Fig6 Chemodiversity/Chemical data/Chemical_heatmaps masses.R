#'--- 
#'Title: Exp3 Beads B.subtilis + syncom --- Chemical
#'Author: Carlos N. Lozano Andrade
#'Aim: 
#'Date: OCt, 22, 2020 
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
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq")
library(vegan)
library(MASS)

dat.Heatmap.d1 <- 
  read.csv('Chemical data/2020.10.22 Chemical_heatmap masses AVERAGE2.csv', header = T, 
           sep = ';', 
           dec =',')
NameVector <- read.csv('Chemical data/NameVector.csv', sep=';', dec = ',')

rownames(dat.Heatmap.d1) <- NameVector$Mass

top50 <- read.csv('Heatmap/Syncom_day4_top50.csv', sep=';', dec = ',')


xx <- NameVector[NameVector$Mass %in% top50$Peaks, ]

df[df$var1 %in% c('value1', 'value2', 'value3'), ]




pheatmap(dat.Heatmap.d1, scale = "row",
         cluster_cols =  F, 
         cluster_rows =  F,
         cellwidth = 12,
         cellheight = 5, 
         height =  8,
         border_color = 'grey')


heatmap.2(dat.Heatmap.d1)



######## Based on Marios selection ############### 


dat.Heatmap.d2 <- 
  read.csv('Chemical data/2020.10.22 Chemical_heatmap masses AVERAGE_clean.csv', header = T, 
           sep = ';', 
           dec =',')
NameVector2 <- read.csv('Chemical data/NameVector2.csv', sep=';', dec = ',')
rownames(dat.Heatmap.d2) <- NameVector2$Mass

pheatmap(dat.Heatmap.d2, scale = "row",
         color = brewer.pal(n = 8, name = "RdYlBu" ),
         cluster_cols =  F, 
         cluster_rows =  T,
         cellwidth = 40,
         cellheight = 10, 
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean", 
         border_color = 'grey', 
         cutree_rows = 5,
         fontsize_row = 12,
         fontsize_col  = 18,
         angle_col = 45)

png ('peaks_heatmap.png', height = 25, width = 45, res = 300, units = 'cm')





data_subset_norm <- t(apply(data_subset, 1, cal_z_score))

my_hclust_gene <- hclust(dist(dat.Heatmap.d2), method = "complete")
plot(my_hclust_gene)
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 3)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))
set.seed(1984)
my_random <- as.factor(sample(x = 1:2, size = nrow(my_gene_col), replace = TRUE))
my_gene_col$random <- my_random
my_sample_col <- data.frame(sample = rep(c("tumour", "normal"), c(4,2)))
row.names(my_sample_col) <- colnames(data_subset)
my_heatmap <- pheatmap(data_subset,
                       annotation_row = my_gene_col,
                       annotation_col = my_sample_col,
                       cutree_rows = 2,
                       cutree_cols = 2)



################# Fulll - all days included ########################### 

dat_full <- read.csv('Chemical data/2020.11.09 Chemical_profile_syncom_long.csv', 
                     sep=';', dec =',')
          ######## Checking the peaks significance ################ 

dat_full$Treatment <- as.factor(dat_full$Treatment)
dat_full$Day <- as.factor(dat_full$Day)
dat_full_sum <- dat_full %>% 
                group_by(Day, Treatment) %>%
                      summarise(
                        across(starts_with('Peak'), list(mean = mean, sd = sd), .names = "{.col}.{.fn}"),
                        na.rm = T
                        ) 


dat_full_sum <- dat_full %>% 
  group_by(Day, Treatment) %>%
  summarise(
    across(starts_with('Peak'), ~mean(.x, na.rm = TRUE))
  ) 




plots <- lapply(colnames(dat_full_sum)[3:length(colnames(dat_full_sum))], function(nm){
  ggplot(dat_full_sum) +
    geom_line(size=2, aes_string(x =colnames(dat_full_sum)[1], colour = colnames(dat_full_sum)[2],
                          nm)) + theme_bw() +
    theme(
      legend.title=element_blank(),axis.text=element_markdown(size=22),
      axis.title=element_text(size=22,face="bold"),
      legend.text = element_markdown(size = 22),
      strip.text.x = element_text(size=22, face="bold")) + 
    scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                 "#87C55F" , "#9EB9F3" )) + 
    scale_x_continuous(breaks = c(1,3, 6, 9)) 
  
})


lapply(names(plots), 
       function(x)ggsave(filename=paste(x,".jpeg",sep=""), plot=plots[[x]]))


for(i in col) {
  print(i)
  
  g<-data_fact %>% group_by(cluster, !!ensym(i)) %>% summarise(total.count=n()) %>% 
    ggplot(., aes(x=cluster, y=total.count, fill=!!ensym(i))) + 
    geom_bar(position  = 'dodge', stat='identity') + 
    geom_text(aes(label=total.count), position = position_dodge(width=0.9), vjust=-0.2) +
    labs(title=i)
  print(g)
}


apply(dat_full[, 4:69], 2, mean, na.rm=T)


apply(dat_full[, 4:69], 2, function(x) tapply(x, dat_full$Day, sum))






#### NMDS ############## 

xfactor <- paste((as.factor(dat_full$Treatment)), (as.factor(dat_full$Day)), sep='_')



peaks.dist <- dist(scale(dat_full[,4:69]))

peaks.cmdscale <- cmdscale(peaks.dist,eig = TRUE,k=2)

# have the class annoation be the coloring guide 
ggplot(as.data.frame(peaks.cmdscale$points), aes(x=peaks.cmdscale$points[,1] , y=peaks.cmdscale$points[,2], 
                                                 col = as.factor(xfactor))) +
  geom_point() +
  labs(title = "classic MDS Plot")




# non-metric MDS
# find the best distance measure
rank.totus <- rankindex(as.factor(xfactor), dat_full[,4:69], indices = c("bray", "euclid", "manhattan", "horn"), method = "spearman")

wines_dist = as.matrix((vegdist(dat_full[,4:69], "manhattan")))


NMDS=isoMDS(wines_dist) # from mass package
NMDS = metaMDS(wines_dist, k=2,trymax=100) # from vegan package
stressplot(NMDS)


#build a data frame with NMDS coordinates and metadata
MDS1 = NMDS$points[,1]
MDS2 = NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, class=as.factor(xfactor))

ggplot(NMDS, aes(x=MDS1, y=MDS2, colour=NMDS$class)) +
  geom_point() +
  labs(title = "NMDS Plot")



# Ordination: Bray_PCoA

nmds = metaMDS(dat_full[,4:69], distance = "bray")
nmds
nmd_factor = envfit(nmds, dat_full[,1:2], permutations = 999, na.rm = TRUE)

plot(nmds)
plot(nmd_factor)


data.scores = as.data.frame(scores(nmds)$sites)
data.scores$Treatment <- dat_full$Treatment
data.scores$Day <- dat_full$Day

en_coord_cont = as.data.frame(scores(nmd_factor, "vectors")) * ordiArrowMul(nmd_factor)
en_coord_cat = as.data.frame(scores(nmd_factor, "factors")) * ordiArrowMul(nmd_factor)



gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(shape = as.factor(Day), colour = Treatment), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                 "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                 "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                 "#968EB9", "#BDDCBC", "#DE9D8B",
                                 "#C0C0C0")) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Treatment", 
       shape = 'Days')



ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(colour = as.factor(Day),  shape = Treatment), size = 3, alpha = 0.5) + 
  scale_colour_manual(values = c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                                 "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                                 "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                                 "#968EB9", "#BDDCBC", "#DE9D8B",
                                 "#C0C0C0"))  + 
  geom_segment(aes(x = 0, y = 0, xend = data.scores$NMDS1, yend = data.scores$NMDS2), 
               data = en_coord_cont, size =1, alpha = 0.5, colour = "grey30") +
  geom_point(data = en_coord_cat, aes(x = NMDS1, y = NMDS2), 
             shape = "diamond", size = 4, alpha = 0.6, colour = "navy") +
  geom_text(data = en_coord_cat, aes(x = NMDS1, y = NMDS2+0.04), 
            label = row.names(en_coord_cat), colour = "navy", fontface = "bold") + 
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Treatment", 
       shape = 'Days')


############################### New NMDS #################################### 


set.seed(19880604)
nmds <- metaMDS(as.matrix(dat_full[,4:69], distance = "bray"))
nmds
plot(nmds)
data.scores = as.data.frame(scores(nmds)$sites)



data.scores = as.data.frame(scores(nmds)$sites)
data.scores$Treatment <- dat_full$Treatment
data.scores$Day <- dat_full$Day



#add columns to data frame 
data.scores$Time = as.factor(dat_full$Day)
data.scores$Treatment = dat_full$Treatment


ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Time, colour = Treatment)) + 
  
  scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                               "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                               "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                               "#968EB9", "#BDDCBC", "#DE9D8B",
                               "#C0C0C0")) + theme_poster() 



+
  geom_segment(data=data.frame(fit_val),
               aes(x=0,y=0,xend=NMDS1, yend=NMDS2),
               arrow=arrow(length=unit(0.2,"cm")),
               color='black',alpha=2/3) +
  geom_label_repel(data=data.frame(fit_val), aes(NMDS1, NMDS2, label=rownames(fit_val)),
                   color='black',alpha=1,
                   segment.color = 'grey35',
                   point.padding = unit(0.3,"lines"),
                   size = 3)


ggsave('Chemistry_data.png',  width=24, height=14, units = "cm", dpi = 600)




######## Removing surfactin ############################ 

vecSurf <- c(9, 15, 16, 17, 53, 58)
data_clean <- dat_full %>% select(!c(15, 16, 17, 53, 58))


set.seed(19880604)
nmdsC <- metaMDS(as.matrix(log(data_clean[,4:6]), distance = "bray"))
nmdsC
plot(nmdsC)
data.scoresC = as.data.frame(scores(nmdsC)$sites)



data.scoresC = as.data.frame(scores(nmdsC)$sites)
data.scoresC$Treatment <- data_clean$Treatment
data.scoresC$Day <- data_clean$Day

fit <- envfit(nmdsC, data_clean[,4:63] ,permutations = 999)
fit_val <- scores(fit, display = c("vectors"))
fit_val <- fit_val*vegan::ordiArrowMul(fit_val, fill = 1.5) 
fit_val <-cbind(fit_val, fit$vectors$pvals)


## fit output:
## p value, Dim1 and Dim2 are vectors added to NMDS plot
fit$vectors


ggplot(data.scoresC, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = Day, colour = Treatment)) + 
  
  scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                               "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                               "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                               "#968EB9", "#BDDCBC", "#DE9D8B",
                               "#C0C0C0")) + theme_poster() +
  geom_segment(data=data.frame(fit_val) %>% filter(V3<0.05),
               aes(x=0,y=0,xend=NMDS1, yend=NMDS2),
               arrow=arrow(length=unit(0.2,"cm")),
               color='black',alpha=2/3) +
  geom_label_repel(data= data.frame(fit_val) %>% filter(V3<0.05),
                   aes(NMDS1, NMDS2, label=rownames(data.frame(fit_val) %>% filter(V3<0.05))),
                   color='black',alpha=1,
                   segment.color = 'grey35',
                   point.padding = unit(0.3,"lines"),
                   size = 3)










fit <- envfit(nmds, dat_full[,4:69] ,permutations = 999)
fit_val <- scores(fit, display = c("vectors"))
fit_val <- fit_val*vegan::ordiArrowMul(fit_val, fill = 1.5)


## fit output:
## p value, Dim1 and Dim2 are vectors added to NMDS plot
fit$vectors






########## PCA##########
library("factoextra")
res.pca <- prcomp(dat_full[,4:69], scale = TRUE)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             geom = c("point", "text"), 
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
)



groups <- as.factor(dat_full$Treatment)
days <- as.factor(dat_full$Day)
fviz_pca_ind(res.pca,
             #fill.ind = days,
             col.ind = groups, # color by groups
             palette = c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                         "#87C55F" , "#9EB9F3"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)   +
  theme(legend.position = "right", 
        legend.title=element_blank(),
        axis.text=element_text(size=22),
        axis.title=element_text(size=22,face="bold"),
        legend.text = element_markdown(size =15),
        strip.text.x = element_text(size=22, face="bold"),axis.ticks.y = element_blank())



ggsave('FOr_poster3.png',  width=18, height=20, units = "cm", dpi = 300) 


daysANOVA <- lapply(split(dat_full, as.factor(dat_full$Day)), 
               function(i){
                    anova(lm(Peak1 ~ Treatment, data = i)) 
                    #agricolae::HSD.test(Peak1 ~ Treatment, data = i)
                 
                 
                 
                        }
               )




daysKrus <- lapply(split(dat_full, as.factor(dat_full$Day)), 
                    function(i){
                     kruskal.test(Peak1 ~ Treatment, data = i) 
                      #agricolae::HSD.test(Peak1 ~ Treatment, data = i)
                    } 
)


daysKb <- lapply(split(dat_full, as.factor(dat_full$Day)), 
                   function(i){
                     boxplot(Peak1 ~ Treatment, data = i) 
                     #agricolae::HSD.test(Peak1 ~ Treatment, data = i)
                   } 
)





dat_list <- split(dat_full, dat_full$Day)
dat_list


anova_results <-  purrr::map(dat_full, ~anova(lm(.x ~ dat_full$Treatment)))




dat_nested <- dat_full %>%
  group_by(Day)  

 xx <- purrr::map(dat_nested[, 4:69], ~anova(lm(.x ~ dat_nested$Treatment)))


anova_results <-
  dat_full %>%
  dplyr::group_by(Day) %>%
  dplyr::do(m1 = summary(aov(Peak1 ~ Treatment, data = .)),
            tuk = HSD.test(aov(Peak1 ~ C, data = .), 'Treatment', console=TRUE ))

dat_full$Day <- as.factor(dat_full$Day)



uniq <- unique(unlist(dat_full$Day))
for (i in 1:length(uniq)){
  data_1 <- subset(data, date == uniq[i])
  #your desired function
}

uniq <- unique(unlist(dat_full$Day))


for ( i in 1:length(uniq) ){
  
  for (j in 4:ncol(dat_full)){
       
       avz <- broom::tidy(aov(dat_full[,j] ~ dat_full$Treatment, 
                              data = subset(dat_full, Day == uniq[i]) ))
       print(avz)
      }
} 
  


#### Day1 ##########

dat_day1 <- dat_full %>% filter (Day == 1) 

res.pca <- prcomp(dat_day1[,4:69], scale = FALSE )


groups <- as.factor(dat_day1$Treatment)
days <- as.factor(dat_full$Day)
fviz_pca_ind(res.pca,
             fill.ind = groups,
             col.ind = groups, # color by groups
             palette = c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                         "#87C55F" , "#9EB9F3"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)

fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)




fviz_pca_ind(res.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(res.pca, repel = TRUE)





apply(mtcars, 2, function(x) tapply(x, mtcars$cyl, mean))






list <- list()
for(i in 4:ncol(dat_day1)){
  
  column <- names(dat_day1[i])
  #tidy will summarise and return neat format
  avz <- broom::tidy(aov(dat_day1[,i] ~ Treatment, data = dat_day1))
 # mult_comp <-HSD.test((aov(dat_day1[,i] ~ Treatment, data = dat_day1)), 'Treatment')$groups
  # Add this condition if you only want aov with P < 0.05 printed
  if(avz$p.value[1] < 0.05) {
    
    #list[i] <- avz$p.value[1]
   print(avz)
   # print(mult_comp)
  }
}

df1 <- do.call(rbind(list), args = F)

p.vals <- sapply(df, function(x){
  
})
df <- dat_day1.1
p.vals <- apply(df, 2, function(x){
  broom::tidy(aov(x ~ Treatment, data = df))
})


dat_day3 <- dat_full %>% filter (Day == 3) 
res.pca <- prcomp(dat_day3[,4:69], scale = FALSE )


groups <- as.factor(dat_day3$Treatment)
days <- as.factor(dat_day3$Day)
fviz_pca_ind(res.pca,
             fill.ind = groups,
             col.ind = groups, # color by groups
             palette = c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                         "#87C55F" , "#9EB9F3"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)

fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)




fviz_pca_ind(res.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(res.pca, repel = TRUE)











for(i in 4:ncol(dat_day3)){
  
  column <- names(dat_day3[i])
  #tidy will summarise and return neat format
  avz <- broom::tidy(aov(dat_day3[,i] ~ Treatment, data = dat_day3))
  mult_comp <-HSD.test((aov(dat_day3[,i] ~ Treatment, data = dat_day3)), 'Treatment')$groups
  # Add this condition if you only want aov with P < 0.05 printed
  if(avz$p.value[1] < 0.05) {
    
    print(column)
    # print(avz)
    print(mult_comp)
    
  }
}



dat_day6 <- dat_full %>% filter (Day == 6) 

res.pca <- prcomp(dat_day6[,4:69], scale = FALSE )


groups <- as.factor(dat_day6$Treatment)
days <- as.factor(dat_day6$Day)
fviz_pca_ind(res.pca,
             fill.ind = groups,
             col.ind = groups, # color by groups
             palette = c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                         "#87C55F" , "#9EB9F3"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)

fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)




fviz_pca_ind(res.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(res.pca, repel = TRUE)
for(i in 4:ncol(dat_day6)){
  
  column <- names(dat_day6[i])
  #tidy will summarise and return neat format
  avz <- broom::tidy(aov(dat_day6[,i] ~ Treatment, data = dat_day6))
  #mult_comp <-HSD.test((aov(dat_day6[,i] ~ Treatment, data = dat_day6)), 'Treatment')$groups
  # Add this condition if you only want aov with P < 0.05 printed
  if(avz$p.value[1] < 0.05) {
    
    print(column)
    # print(avz)
    #print(mult_comp)
    
  } 
}


dat_day9 <- dat_full %>% filter (Day == 9) 

res.pca <- prcomp(dat_day9[,4:69], scale = FALSE )


groups <- as.factor(dat_day9$Treatment)
days <- as.factor(dat_day6$Day)
fviz_pca_ind(res.pca,
             fill.ind = groups,
             col.ind = groups, # color by groups
             palette = c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                         "#87C55F" , "#9EB9F3"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)

fviz_pca_var(res.pca, col.var="contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping
)




fviz_pca_ind(res.pca, col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)
fviz_pca_biplot(res.pca, repel = TRUE)

for(i in 4:ncol(dat_day9)){
  
  column <- names(dat_day9[i])
  #tidy will summarise and return neat format
  avz <- broom::tidy(anova(lm(dat_day9[,i] ~ Treatment, data = dat_day9)))
  #mult_comp <-HSD.test((aov(dat_day6[,i] ~ Treatment, data = dat_day6)), 'Treatment')$groups
  # Add this condition if you only want aov with P < 0.05 printed
  if(avz$p.value[1] < 0.05) {
    
    print(column)
    # print(avz)
    #print(mult_comp)
    
  } 
}




dat_day1.1 <- dat_full %>% filter (Day == 1)  %>%  select(-c('Day', 'Replicate'))



aov.modelsD1 <- lapply(setdiff(names(dat_day1.1), "Treatment"), function(s) {
  anova(lm(as.formula(paste(s, " ~ Treatment")),dat_day1.1))
})


dat_day3.1 <- dat_full %>% filter (Day == 3)  %>%  select(-c('Day', 'Replicate'))



aov.modelsD3.1 <- lapply(setdiff(names(dat_day3.1), "Treatment"), function(s) {
  anova(lm(as.formula(paste(s, " ~ Treatment")),dat_day3.1))
})




dat_day6.1 <- dat_full %>% filter (Day == 6)  %>%  select(-c('Day', 'Replicate'))


aov.modelsD6.1 <- lapply(setdiff(names(dat_day6.1), "Treatment"), function(s) {
  anova(lm(as.formula(paste(s, " ~ Treatment")),dat_day6.1))
})




dat_day9.1 <- dat_full %>% filter (Day == 9)  %>%  select(-c('Day', 'Replicate'))


aov.modelsD9.1 <- lapply(setdiff(names(dat_day9.1), "Treatment"), function(s) {
  anova(lm(as.formula(paste(s, " ~ Treatment")),dat_day9.1))
})




hsD1 <- lapply(setdiff(names(dat_day1), "Treatment"), function(s) {
  agricolae::HSD.test((lm(as.formula(paste(s, " ~ Treatment")),dat_day1)), 'Treatment' )
    })



############################ Heatmpas by day with data full##################### 


dat.Heatmap.d1 <- 
  read.csv('Chemical data/heatmap_day1_full.csv', header = T, 
           sep = ';', 
           dec =',')
NameVectorF <- read.csv('Chemical data/NameVector_Full.csv', sep=';', dec = ',')
rownames(dat.Heatmap.d1) <- NameVectorF$Mass

pheatmap(dat.Heatmap.d1, scale = "row",
         color = brewer.pal(n = 8, name = "RdYlBu" ),
         cluster_cols =  F, 
         cluster_rows =  T,
         cellwidth = 20,
         cellheight = 5, 
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean", 
         border_color = 'grey', 
         cutree_rows = 5,
         fontsize_row = 10,
         fontsize_col  = 18,
         angle_col = 45)

png('peaks_heatmap_day1.png', height = 25, width = 45, res = 300, units = 'cm')



################## Surfactin dynamic################ 

surf_dynamic <- read.csv('Chemical data/Surfactin_dynamic.csv',  header = T, 
                         sep = ';', 
                         dec =',')
legend_labels <- c( "*ppsC*","*WT*")

ggplot(surf_dynamic, aes(x=Day, y=Concentration*5, colour=Strain)) + 
  geom_point() + geom_smooth()  +  ylim(c(0, 50*5)) +
  labs ( y = 'Concentration (µg/mL)', x = 'Time (h)') + theme_bw() +
  theme(
    legend.title=element_blank(),axis.text=element_markdown(size=22),
    axis.title=element_text(size=22,face="bold"),
    legend.text = element_markdown(size = 22),
    strip.text.x = element_text(size=22, face="bold")) + 
  scale_colour_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                               "#87C55F" , "#9EB9F3" ), labels = legend_labels) + 
  scale_x_continuous(breaks = c(1,3, 6, 9))

ggsave ('surfactin_plot.png', height = 15, width = 25, dpi = 300, units = 'cm')
