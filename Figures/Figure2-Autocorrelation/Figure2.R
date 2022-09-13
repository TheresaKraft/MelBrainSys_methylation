setwd(paste0(figureDirectory, "Figure2-Autocorrelation" ) )

library(ggplot2)
load(file = "autocorrelationData.RData")

autocorrelationQuantilesData <- autocorrelationData[[1]]
autocorrelationQuantilesPermut <- autocorrelationData[[2]]

autocorrelationQuantilesData$lag <- autocorrelationQuantilesData$lag -1 


pdf(file ="Figure2.pdf", width = 7 , height = 5)
ggplot(data = autocorrelationQuantilesData[-1,], mapping=aes(x=lag)) +
  geom_ribbon(aes(ymin = quant_0, ymax = quant_100, fill = "0%-100%"), alpha = .5) +
  geom_ribbon(aes(ymin = quant_10, ymax = quant_90, fill = "10%-90%"), alpha = .5) +
  geom_ribbon(aes(ymin = quant_25, ymax = quant_75, fill = "25%-75%"), alpha = .5) +
  geom_point(mapping = aes(x=lag, y = quant_50, colour = "Original data")) +   geom_line(mapping = aes(x=lag, y = quant_50, colour = "Original data")) + 
  geom_line(mapping = aes(x = lag, y = unlist(autocorrelationQuantilesPermut[-1,]$quant_50), colour = "Permuted data")) + 
  geom_point(mapping = aes(x=lag, y = unlist(autocorrelationQuantilesPermut[-1,]$quant_50), colour = "Permuted data")) + 
  scale_x_continuous( breaks = c(1,25,50,75,100), lim = c(0,100)) +
  geom_ribbon(aes(ymin = unlist(autocorrelationQuantilesPermut[-1,]$quant_0), ymax = unlist(autocorrelationQuantilesPermut[-1,]$quant_100), fill = "0%-100%" ), alpha = .5) +
  scale_fill_manual(values=c("aquamarine4","cornflowerblue", "chocolate"), name="Quantiles") +
  scale_color_manual(values=c("red", "black")) +  ylab("ACF") + labs(x = "Positional lag of neighboring CpGs", y = "Autocorrelation", fill = "Quantiles", colour = "Median")
dev.off()
