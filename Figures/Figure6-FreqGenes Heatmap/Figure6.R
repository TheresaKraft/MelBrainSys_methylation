setwd(paste0(figureDirectory, "Figure6-FreqGenes Heatmap"))

library(ggplot2)
library(pheatmap)
library(grid)
library(stringr)


load("mostFreqAnnotationRowList.RData")


mostFreqCpGs <-mostFreqAnnotationRowList[[1]]
annotationRows <- mostFreqAnnotationRowList[[2]]


ann_colorsN = list(  
  cytokine = c("TRUE" = "darkorange", "FALSE" = "lightgrey"),
  Pi3kAkt = c("TRUE" = "deeppink", "FALSE" = "lightgrey"),
  ECM = c("TRUE" = "chartreuse", "FALSE" = "lightgrey"),
  MAPK = c("TRUE" = "aquamarine2", "FALSE" = "lightgrey"),
  Region = c(body = "bisque3", promoter = "blueviolet", enhancer = "darkolivegreen3")
)


htmp <- pheatmap(mostFreqCpGs, color = colorRampPalette(c('cornflowerblue','white','brown3'))(length(seq(-6,6, by = 0.01))), breaks = seq(-6,6, by = 0.01), 
                 clustering_distance_rows = "manhattan", clustering_distance_cols = "manhattan", legend_breaks = c(-4, 4), 
                 main = "", legend_labels = c("Decreased \n methylation", "Increased \n methylation"),
                 annotation_row  = annotationRows, annotation_color =  ann_colorsN,annotation_legend	=F,  show_rownames=T, treeheight_row = F, fontsize_row = 7)




labelshtmp <- htmp$tree_col$labels[htmp$tree_col$order]
labelshtmp[str_detect(labelshtmp ,"BLun")] <- "blue3"
labelshtmp[str_detect(labelshtmp , "BLym")] <- "forestgreen"
labelshtmp[str_detect(labelshtmp , "BSki")] <- "gold3"
labelshtmp[str_detect(labelshtmp , "BSof")] <- "orchid"
labelshtmp[str_detect(labelshtmp , "BLiv")] <- "deeppink2"

cols= labelshtmp
htmp$gtable$grobs[[4]]$gp=gpar(col=cols)




pdf(file ="Figure6.pdf",  width =7 , height = 8)
htmp
dev.off()
