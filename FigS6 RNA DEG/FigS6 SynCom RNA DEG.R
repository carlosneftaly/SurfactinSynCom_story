
data.plot <-read.csv('up_down_genes.csv', header = T, sep = ';')



data.plot$direction <- factor(data.plot$direction, 
                       levels = c("Up-regulated", "Down-regulated"), 
                       labels = c('Up-regulated',
                                  'Down-regulated')
)
ggplot(data.plot, aes(x=Strain, y=n.DEGs, fill=direction)) + theme_thesis() +
  geom_bar(stat="identity", position="identity") + facet_grid(~Day) +
  theme(legend.position="bottom", legend.title=element_blank()) +
  labs(y='# of DE genes') +  
  scale_y_continuous(labels = c(200, 100, 0, 100), expand = c(0,10))  +
  scale_fill_manual(values=c("#79B4B8", "#990000", "#F6CF71","#DCB0F2", 
                             "#87C55F" , "#9EB9F3" ,"#FE88B1" , "#C9DB74","#FEDDC4",
                             "#8BE0A4", "#B497E7" ,"#D3B484", "#E2F0CB", 
                             "#968EB9", "#BDDCBC", "#DE9D8B",
                             "#C0C0C0"))  

ggsave('DEG_summary_full.pdf', width = 12, height = 10, dpi = 600, units = 'in')

ggsave('DEG_summary_full.png', width = 12, height = 10, dpi = 300, units = 'in')

