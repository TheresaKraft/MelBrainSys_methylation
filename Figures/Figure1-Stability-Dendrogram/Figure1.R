# Set directory to source file loaction
setwd(paste0(figureDirectory, "Figure1-Stability-Dendrogram" ) )

library(RColorBrewer)
library(stringr)
library(pvclust)
library(dendextend)

load(file = "stabilityHclust_wardD2_manhattan.RData")


patientColorManhD2 <- valuesManhD2.pv$hclust$labels[valuesManhD2.pv$hclust$order]
patientColorManhD2[str_detect(patientColorManhD2,"P17")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[1]
patientColorManhD2[str_detect(patientColorManhD2,"P16")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[2]
patientColorManhD2[str_detect(patientColorManhD2,"P09")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[3]
patientColorManhD2[str_detect(patientColorManhD2,"P67")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[4]
patientColorManhD2[str_detect(patientColorManhD2,"P39")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[5]
patientColorManhD2[str_detect(patientColorManhD2,"P28")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[6]
patientColorManhD2[str_detect(patientColorManhD2,"P03")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[7]
patientColorManhD2[str_detect(patientColorManhD2,"P42")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[8]
patientColorManhD2[str_detect(patientColorManhD2,"P18")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[9]
patientColorManhD2[str_detect(patientColorManhD2,"P04")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[10]
patientColorManhD2[str_detect(patientColorManhD2,"P08")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[11]
patientColorManhD2[str_detect(patientColorManhD2,"P64")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[12]
patientColorManhD2[str_detect(patientColorManhD2,"P02")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[13]
patientColorManhD2[str_detect(patientColorManhD2,"P06")] <- colorRampPalette(brewer.pal(8, "Oranges"))(14)[14]

tissueColorManhD2 <- valuesManhD2.pv$hclust$labels[valuesManhD2.pv$hclust$order]
tissueColorManhD2[str_detect(tissueColorManhD2,"_B")] <- "darkgrey"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Lym")] <- "forestgreen"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Lun")] <- "blue3"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Liv")] <- "deeppink2"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Ski")] <- "gold3"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Sof")] <- "darkorchid"



pdf(file ="Figure1.pdf", width = 9, height = 5.50 )
par(mar = c(3,3,0,0))
plot(valuesManhD2.pv, hang = -1, cex = 0.9, print.pv = c("au"), main ="", axes = F, labels =  c("P42_B","P42_Lym-1","P42_Lym-2","*P39_B",                                                                                            "*P39_Lun","P28_B","P28_Ski","P18_B","P18_Lun-1","P18_Lun-2","*P17_B","*P17_Lym-1","*P17_Lym-2",
                                                                                                "P06_B-1","P06_B-2","P06_Lym","P04_B","P04_Ski-1","P04_Ski-2","P03_B","P03_Lun","P02_B","P02_Lym",
                                                                                                "*P09_B","*P09_Liv","P08_B","P08_Sof-1","P08_Sof-2","P08_Sof-3","*P67_B-1","*P67_B-2","*P67_Lun-1",
                                                                                                "*P67_Lun-2","P64_B","P64_Lun","P16_B","P16_Lun"),
     print.num = F)
colored_bars(colors = cbind(patientColorManhD2, tissueColorManhD2),rowLabels = c("Patient","Tissue"))
dev.off()


