setwd(paste0(figureDirectory, "Figure7-expression"))

library(ggplot2)
library(grid)


load("exprMethylConnection.RData")
load("ExpressionComparisonPlot.RData")

plotExpressionMethylation <- ggplot(exprMethylConnection, aes (x=reorder(Gene,medianDiff), y = value, fill = data)) + 
  geom_hline(yintercept = 0, col = "black") + geom_boxplot() +
  labs(x = "", y = "log2-ratio( intra / extra )", fill = "Data") + 
  scale_color_manual(values = c("bisque3", "red"), guide = "none") + scale_fill_manual(values = c("peachpuff1", "salmon2")) + 
  geom_vline(xintercept = c(8.5,11.5), col = "black") +  scale_y_continuous(breaks = seq(-5, 8, by = 1)) + 
  theme(legend.position = c(0.945,0.88), legend.background = element_rect(size=0.5, linetype="solid",
                                                                         colour ="black"))  + labs(fill = "Measurements")


# gene TNXB

t.test(ExpressionComparisonPlot$value[ExpressionComparisonPlot$gene == "TNXB" &ExpressionComparisonPlot$region == "intracranial" ], 
       ExpressionComparisonPlot$value[ExpressionComparisonPlot$gene == "TNXB" &ExpressionComparisonPlot$region == "extracranial" ],
       paired = T)


plotTNXB <- ggplot(ExpressionComparisonPlot[ExpressionComparisonPlot$gene == "TNXB",], aes(x = region, y = value, fill = region )) + geom_boxplot() + 
  labs(x="TNXB", y = "log2 expression value")  + scale_fill_manual(values = c("white","darkgrey")) + 
  theme(legend.position = "none") + scale_y_continuous(limits = c(0,10), breaks = seq(0,10,1)) + 
  geom_text(aes(x = 1.5, y = 9), label = "p = 0.031" )

# gene JAK3

t.test(ExpressionComparisonPlot$value[ExpressionComparisonPlot$gene == "JAK3" &ExpressionComparisonPlot$region == "intracranial" ], 
       ExpressionComparisonPlot$value[ExpressionComparisonPlot$gene == "JAK3" &ExpressionComparisonPlot$region == "extracranial" ],
       paired = T)


plotJAK3 <- ggplot(ExpressionComparisonPlot[ExpressionComparisonPlot$gene == "JAK3",], aes(x = region, y = value, fill = region )) + geom_boxplot() + 
  labs(x="JAK3", y = "log2 expression value")   + scale_fill_manual(values = c("white","darkgrey")) + 
  theme(legend.position = "none") + scale_y_continuous(limits = c(0,10), breaks = seq(0,10,1)) + 
  geom_text(aes(x = 1.5, y = 9), label = "p = 0.048" )

# gene MECOM

t.test(ExpressionComparisonPlot$value[ExpressionComparisonPlot$gene == "MECOM" &ExpressionComparisonPlot$region == "intracranial" ], 
       ExpressionComparisonPlot$value[ExpressionComparisonPlot$gene == "MECOM" &ExpressionComparisonPlot$region == "extracranial" ],
       paired = T)


plotMECOM <- ggplot(ExpressionComparisonPlot[ExpressionComparisonPlot$gene == "MECOM",], aes(x = region, y = value, fill = region )) + geom_boxplot() + 
  labs(x="MECOM", y = "log2 expression value")  + geom_vline(xintercept = c(10.5,15.5,16.5,19.5)) + scale_fill_manual(values = c("white","darkgrey")) + 
  theme(legend.position = "none") + scale_y_continuous(limits = c(0,10), breaks = seq(0,10,1)) + 
  geom_text(aes(x = 1.5, y = 9), label = "p = 0.033" )




lower <- plot_grid(plotTNXB, plotMECOM, plotJAK3, nrow = 1,  label_size = 15,  labels = c('B', 'C','D')) 
finalPlot <- plot_grid(plotExpressionMethylation, lower,  nrow = 2, rel_heights = c(60,40), label_size = 15, labels = c('A',''))



pdf(file ="Figure7.pdf",   width = 11 , height = 7)
finalPlot
dev.off()
