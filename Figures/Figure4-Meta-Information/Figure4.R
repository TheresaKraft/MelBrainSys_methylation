setwd(paste0(figureDirectory, "Figure4-Meta-Information"))

library(ggplot2)
library(cowplot)

load("metaInfoDF.RData")


metaInfoDF$alive <- ""
metaInfoDF$alive[metaInfoDF$Patient == "P39-Blun"] <- "*"


plotLegend <- ggplot(metaInfoDF, aes(shape=TherapyTissue,y=survivalBM , x=numDifCpG,color = negScore ) ) +
  geom_point(cex = 4) +  labs(color='Methylation \n states', shape = "Metastasis treated") + theme(legend.direction = "vertical",  legend.box = "vertical") +
  scale_color_gradient2(low = "cornflowerblue", mid = "azure4", high = "brown3", breaks = c(-0.5,0.4), 
                          labels = c("More \n decreased \n CpGs","More \n increased \n CpGs"))  + theme(legend.direction = "vertical", 
                                                                                                        legend.box = "vertical")


allLegend <- get_legend(plotLegend)

plot1 <- ggplot(metaInfoDF, aes(shape=TherapyTissue,y=timeBetween , x=numDifCpG , color = negScore) ) +
  labs(color='Number of \n altered CpGs', shape = "Therapy", 
       x ="Number of differentially methylated CpGs \n intra- vs. extracranial", 
       y = "Time between metastases \n (months)") +
  theme(legend.title = element_blank(), legend.position = "none") + 
  geom_point(cex = 3) + scale_color_gradient2(low = "cornflowerblue", mid = "azure4", high = "brown3", breaks = c(-0.5,0.4), 
                                              labels = c("More \n decreased \n CpGs","More \n increased \n CpGs"))  

  
  
plot2 <- ggplot(metaInfoDF, aes(shape=TherapyTissue,y= therapyOption, x=numDifCpG , color = negScore) ) +
  labs(color='Number of \n altered CpGs', shape = "Therapy", 
       x ="Number of differentially methylated CpGs \n intra- vs. extracranial", 
       y = "") +
  theme(legend.title = element_blank(), legend.position = "none") + 
  geom_point(cex = 3) +   scale_color_gradient2(low = "cornflowerblue", mid = "azure4", high = "brown3", breaks = c(-0.5,0.4), 
                                                labels = c("More \n decreased \n CpGs","More \n increased \n CpGs")) 



plot3 <- ggplot(metaInfoDF, aes(shape=TherapyTissue,y=survivalBM , x=numDifCpG , color = negScore, label = alive) ) +
  labs(color='Number of \n altered CpGs', shape = "Therapy", 
       x ="Number of differentially methylated CpGs \n intra- vs. extracranial", 
       y = "\n Survival from intracranial metastasis \n (months) ") +
  theme(legend.title = element_blank(), legend.position = "none") + 
  geom_point(cex = 3) +   scale_color_gradient2(low = "cornflowerblue", mid = "azure4", high = "brown3", breaks = c(-0.5,0.4), 
                                                labels = c("More \n decreased \n CpGs","More \n increased \n CpGs")) + geom_text(hjust=0, vjust=0, color = "black")



pdf(file ="Figure4.pdf", width = 12 , height = 4)
plot_grid(plot1, plot2, plot3, allLegend, align = "h", nrow = 1, rel_widths = c(28/100,30/100,28/100,14/100), labels = c('A', 'B','C',''), label_size = 18)
dev.off()
  
  
