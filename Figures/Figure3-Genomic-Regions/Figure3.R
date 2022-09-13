setwd( paste0(figureDirectory, "Figure3-Genomic-Regions") )

library(ggplot2)
library(cowplot)


load("barplotOverviewDF.RData")
load("enhancerGenePlotDf.RData")
load("pathPlotDfAllFreq.RData")

colors <- c(rep("blue3",10),rep("forestgreen",7),rep("gold3",3),rep("darkorchid",3 ),rep("deeppink2",1))


plot1 <- ggplot(barplotOverviewDF, aes(color=Methylation,y=value, x=reorder(PatientID,-ranking)), ylim = c(0,4)) + 
  geom_bar(position = "identity",stat="identity", fill="white",) +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill="white", size=.5, linetype="solid",color="black"),
        legend.title = element_blank(), legend.direction = "horizontal") +
  scale_y_continuous(breaks = c(-2e+05,-1e+05,0,1e+05,2e+05), labels = c("2e+05","1e+05","0","1e+05","2e+05"), limits = c(-2.1e+05,2.1e+05))+
  ggtitle("") + 
  ylab("Number of CpGs \n Decreased  \t \t Increased") +
  geom_vline(xintercept = c(10.5,17.5,20.5,23.5)) +  geom_hline(yintercept = 0)   +  scale_color_manual(values= c("cornflowerblue","brown3") ) 

enhancerGenePlotDf <- enhancerGenePlotDf[enhancerGenePlotDf$association != "Others",]

plot2 <- ggplot(enhancerGenePlotDf, aes(fill=association,y=value, x=reorder(PatientID,-ranking)), ylim = c(0,4)) + 
  geom_bar(width = 0.7, stat="identity", position = position_dodge(width = 0.7), colour="white") + 
  theme(axis.text.x =element_blank(),
        axis.title.x=element_blank(),legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill="white", size=.5, linetype="solid",color="black"),
        legend.title = element_blank(),legend.direction = "horizontal") +
  ggtitle("") +  ylab(" \n % affected CpGs \n Decreased \t \t Increased") + 
  geom_text(aes(label=significancePos),stat='identity',position=position_dodge(0.68),vjust=0, size = 3, colour = "red",fontface = "bold")  +
  geom_text(aes(label=significanceNeg),stat='identity',position=position_dodge(0.68),vjust=1, size = 3, colour = "red",fontface = "bold")  +
  scale_y_continuous(breaks = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4), labels = c("40","30","20","10","0","10","20","30","40"), limits = c(-0.30,0.30))+
  geom_vline(xintercept = c(10.5,17.5,20.5,23.5)) +  geom_hline(yintercept = 0) + scale_fill_manual(values= c("bisque3","blueviolet","darkolivegreen3") )


plot3 <- ggplot(pathPlotDfAllFreq, aes(fill=association,y=value, x=reorder(PatientID,-ranking)), ylim = c(0,4)) + 
  geom_bar(width = 0.7, stat="identity", position = position_dodge(width = 0.7), colour="white") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,color = colors),
        axis.title.x=element_blank(),legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill="white", size=.5, linetype="solid",color="black"),
        legend.title = element_blank(),legend.direction = "horizontal") +
  ggtitle("") +  ylab("% affected CpGs \n Decreased  \t \t Increased ") + 
  geom_text(aes(label=significancePos),stat='identity',position=position_dodge(0.68),vjust=0, size = 3, colour = "red")  +
  geom_text(aes(label=significanceNeg),stat='identity',position=position_dodge(0.68),vjust=1, size = 3, colour = "red")  +
  scale_y_continuous(breaks = c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4), labels = c("40","30","20","10","0","10","20","30","40"), limits = c(-0.35,0.3))+
  geom_vline(xintercept = c(10.5,17.5,20.5,23.5)) +  geom_hline(yintercept = 0) + scale_fill_manual(values= c("darkorange","chartreuse","cornflowerblue",
                                                                                                              "deeppink","aquamarine2","deepskyblue4"))

pdf(file ="Figure3.pdf", width = 12, height = 10)
plot_grid(plot1, plot2, plot3, align = "v", nrow = 3, rel_heights = c(31/100, 31/100, 38/100), label_size = 15, labels = c('A', 'B', 'C')) 
dev.off()

