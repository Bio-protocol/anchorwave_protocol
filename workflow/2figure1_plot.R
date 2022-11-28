#Using R to generate a dotplot.
library(ggplot2)
#Transform Coordinates using a function.
changetoM <- function ( position ){
position=position/1000000;
paste(position, "M", sep="")
}
#Read gene position, belong to which chromosome and so on
data =read.table("cds.tab")
#Select all euchromosomes as factor
data = data[which(data$V1 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data = data[which(data$V3 %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")),]
data$V1 = factor(data$V1, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
data$V3 = factor(data$V3, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))
#Using ggplot2 to plot a dotplot and beautify it.
png("E:\\R_Scripts\\figure1.png")
ggplot(data=data, aes(x=V4, y=V2)) +geom_point(size=0.5, aes(color=V5)) +
facet_grid(V1 ~ V3, scales="free",space="free") +labs(x="sorghum", y="maize")+
scale_x_continuous(labels=changetoM) + scale_y_continuous(labels=changetoM) +
theme(axis.line = element_blank(),
panel.background = element_blank(),
panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
axis.text.y = element_text( colour = "black"),
legend.position='none',
axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
dev.off()
