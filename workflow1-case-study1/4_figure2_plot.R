library(ggplot2)
library(svglite)
changetoM <- function ( position ){
position=position/1000000;
paste(position, "M", sep="")}
data =read.table("align1.anchors", header=TRUE)
data = data[which(data$refChr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data = data[which(data$queryChr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data$refChr = factor(data$refChr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
data$queryCh = factor(data$queryChr, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
figure2 <- ggplot(data=data, aes(x=queryStart, y=referenceStart))+
geom_point(size=0.5, aes(color=strand)) +
facet_grid(refChr~queryChr, scales="free", space="free") +
labs(x="sorghum", y="maize")+scale_x_continuous(labels=changetoM) +
scale_y_continuous(labels=changetoM) +
  theme(axis.line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill =NA,color="black", size=0.5, linetype="solid"),
        axis.text.y = element_text( colour = "black"),
        legend.position='none',
        axis.text.x = element_text(angle=300,hjust=0, vjust=0.5, colour = "black") )
png("figure2.png")  
figure2
dev.off()
pdf("figure2.pdf")
figure2
dev.off()
svglite("figure2.svg")
figure2
dev.off()
