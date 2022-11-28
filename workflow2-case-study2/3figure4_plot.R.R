setwd("E:\\plot\\maize")
#Use this piece of R to draw a dotplot.
library(ggplot2)
#Transform Coordinates using a function.
changetoM <- function ( position ){
  position=position/1000000;
  paste(position, "M", sep="")
}
#Read gene position, belong to which chromosome and so on
data =read.table("align.anchors", header=TRUE)
#Select all euchromosomes as factor, modify the chromosomes number manually.
data = data[which(data$refChr %in% c("1", "2", "3", "4",
                                     "5", "6", "7", "8", "9", "10")),]
data = data[which(data$queryChr %in% c("1", "2", "3", "4",
                                       "5", "6", "7", "8", "9", "10")),]
data$refChr = factor(data$refChr, levels=c("1", "2", "3", "4",
                                           "5", "6", "7", "8", "9", "10"))
data$queryCh = factor(data$queryChr, levels=c("1", "2", "3", "4",
                                              "5", "6", "7", "8", "9", "10"))
#Using ggplot2 to plot a dotplot and beautify it.
png("E:\\plot\\maize\\figure4.png")
ggplot(data=data, aes(x=queryStart, y=referenceStart))+
  geom_point(size=0.5, aes(color=strand)) +
  facet_grid(refChr~queryChr, scales="free", space="free") +
  labs(x="Mo17", y="B73")+scale_x_continuous(labels=changetoM) +
  scale_y_continuous(labels=changetoM) +
  theme(axis.line = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill =NA,color="black", size=0.5, linetype="solid"),
        axis.text.y = element_text( size=5,colour = "black"),
        legend.position='none',
        axis.text.x = element_text(angle=300,size=6,hjust=0, vjust=0.5, colour = "black") )
dev.off()
