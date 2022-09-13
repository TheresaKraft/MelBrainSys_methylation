setwd(paste0(figureDirectory, "Figure5-genome-wide-candidates" ) )

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(grid)


load("numCand.RData")
load("compositionPathDF.RData")
load("compositionCancDF.RData")
load("enrichdf2.RData")

numCand <- numCand[numCand$Level == "Gene",]
enriched <- enrichdf2[enrichdf2$qval > 2,]
enriched$perc <- as.numeric(enriched$perc)

compositionPathDF$percentage <- as.numeric(as.character(compositionPathDF$percentage))
compositionPathDF <- na.omit(compositionPathDF[compositionPathDF$percentage > 0 , ])

compositionPathDF$pathway <- as.character(compositionPathDF$pathway)
compositionCancDF$cancerPath <- as.character(compositionCancDF$cancerPath)

compositionPathDF$logQVal <- -log10(compositionPathDF$qvalFisher)
compositionCancDF$logQVal <- -log10(compositionCancDF$qvalFisher)

compositionPathDF$pathway[compositionPathDF$pathway == "erbb"] <- "ErbB"
compositionPathDF$pathway[compositionPathDF$pathway == "mapk"] <- "MAPK"
compositionPathDF$pathway[compositionPathDF$pathway == "pi3k_akt"] <- "Pi3K/Akt"
compositionPathDF$pathway[compositionPathDF$pathway == "mtor"] <- "mTOR"
compositionPathDF$pathway[compositionPathDF$pathway == "adh_jun"] <- "adherens junction"
compositionPathDF$pathway[compositionPathDF$pathway == "wnt"] <- "Wnt"
compositionPathDF$pathway[compositionPathDF$pathway == "tgfb"] <- "TGFB"
compositionPathDF$pathway[compositionPathDF$pathway == "focal_adh"] <- "focal adhesion"
compositionPathDF$pathway[compositionPathDF$pathway == "jak_stat"] <- "JAK-STAT"
compositionPathDF$pathway[compositionPathDF$pathway == "cell_cycle"] <- "cell cycle"
compositionPathDF$pathway[compositionPathDF$pathway == "vegf"] <- "VEGF"
compositionPathDF$pathway[compositionPathDF$pathway == "NER"] <- "DNA repair"

compositionCancDF$cancerPath[compositionCancDF$cancerPath == "CancerCensusGene"] <- "cancer census"
compositionCancDF$cancerPath[compositionCancDF$cancerPath == "Kinase"] <- "kinase"
compositionCancDF$cancerPath[compositionCancDF$cancerPath == "Phosphatase"] <- "phosphatase"
compositionCancDF$cancerPath[compositionCancDF$cancerPath == "Oncogene"] <- "oncogene"
compositionCancDF$cancerPath[compositionCancDF$cancerPath == "TF.Cofactor"] <- "transcription factor"
compositionCancDF$cancerPath[compositionCancDF$cancerPath == "TumorSuppressorGene"] <- "tumor suppressor"
compositionCancDF$cancerPath[compositionCancDF$cancerPath == "SignalingPathwayGene"] <- "signaling pathway"
compositionCancDF$cancerPath[compositionCancDF$cancerPath == "MetabolomePathwayGene"] <- "metabolome pathway"
compositionCancDF$cancerPath[compositionCancDF$cancerPath == "EssentialGene"] <- "essential gene"

compositionCancDF$significance <- factor(compositionCancDF$qvalFisher <0.05, levels = c(TRUE, FALSE), labels = c("X",""))
compositionPathDF$significance <- factor(compositionPathDF$qvalFisher <0.05, levels = c(TRUE, FALSE), labels = c("X",""))

compositionPathDF$pathway <- as.factor(compositionPathDF$pathway)
compositionCancDF$cancerPath <- as.factor(compositionCancDF$cancerPath)


legendPositionPercaffected <-  c(0.8,0.4)
legendPositionGO <- c(0.8,0.08)

yAxisName <- expression(-log[10](q))

plot1 <- ggplot(numCand, aes(Patients, Number, color = Direction)) + geom_point(cex = 3) +  geom_line(aes(Patients,Number, color = Direction, shape = Level)) +
  scale_y_continuous(trans = "log10", name = "Number of gene candidates", breaks = c(1,10,100,1000,10000,10000,100000), labels = c(1,10,100,1000,10000,10000,"100000")) + 
  scale_x_continuous(name = "Number of patients cutoff", breaks = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14), minor_breaks = NULL)  +
  scale_color_manual(values = c("cornflowerblue","brown3")) + theme(legend.position = "none")+ 
  theme(legend.position = c(0.915,0.79), legend.background = element_rect(size=0.5, linetype="solid", 
                                                                     colour ="black"))  + labs(color = "Intracranial \nmethylation")


plot2 <- ggplot(compositionPathDF, aes(x = reorder(pathway,percentage), y =  percentage)) + geom_bar(stat='identity', col = "black", fill = "darkgrey") + 
  scale_y_continuous(name = "% pathway genes", limits = c(0,6)) +
  xlab("") + coord_flip()  + scale_fill_manual(values  = c("red","darkgrey"))+  
  theme(legend.position = "none", legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) 


plot3 <- ggplot(compositionCancDF, aes(x = reorder(cancerPath,percentage), y = percentage, fill = significance )) + geom_bar(stat='identity', col = "black") + 
  scale_y_continuous(name = "% pathway genes", limits = c(0,3)) + 
  xlab("") + coord_flip()  +scale_fill_manual(values  = c("red","darkgrey"))+ 
  theme(legend.position = "none", legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) 

plotLegend <- ggplot(compositionPathDF, aes(x = reorder(pathway,percentage), y = percentage )) + 
  geom_bar(stat='identity', col = "black") + 
  scale_y_continuous(name = yAxisName) + 
  xlab("") + coord_flip() + scale_fill_gradientn(colors = brewer.pal(8,"Reds"),limits=c(0,9),  breaks=c(0,3,6,9),name = "% affected genes") + 
  theme(legend.position = "bottom", legend.background = element_rect(size=0.5, linetype="solid", colour ="black")) 
legendPlots <- get_legend(plotLegend)

plotGO <- ggplot(enriched, aes(x = reorder(Description,qval), y = qval, fill = factor(GOType, levels = c("Function", "Composition"), labels = c("MF","CC")))) +
  geom_bar(stat='identity', col = "black") + 
  scale_y_continuous(name = yAxisName, breaks = c(0,1,2,3,4,5)) + labs(fill="Ontology") +
  xlab("") + coord_flip() + scale_fill_manual(values  = brewer.pal(2,"Dark2")) + 
  theme(legend.position = "bottom", legend.background = element_rect(size=0.5, linetype="solid", 
                                                                      colour ="black"))
legendGO <- get_legend(plotGO)

plot4 <- ggplot(enriched, aes(x = reorder(Description,qval), y = qval, fill = factor(GOType, levels = c("MF", "CC"), 
                                                                                     labels = c("Molecular \n function","Cellular \n composition")))) +
  geom_bar(stat='identity', col = "black") + 
  scale_y_continuous(name = yAxisName, breaks = c(0,1,2,3,4,5), limits = c(0,5)) + labs(fill="Ontology") +
  xlab("") + coord_flip() + scale_fill_manual(values  = brewer.pal(2,"Dark2"))  + 
  theme(legend.position = "none", legend.background = element_rect(size=0.5, linetype="solid", 
                                                                     colour ="black"))


middlePlot <- plot_grid(plot2, plot3, nrow = 1,  label_size = 15,  labels = c('B', 'C')) 
leftPlot <- plot_grid(plot1, middlePlot, nrow = 2,  label_size = 15,  labels = c('A', ' '))
finalPlot <- plot_grid(leftPlot, plot4,  nrow = 1,rel_widths = c((60/100),(40/100)), label_size = 15, labels = c('','D'))




pdf(file ="Figure5.pdf",   width = 11 , height = 6)
finalPlot
dev.off()
