#!/usr/bin/env Rscript
library(dplyr)
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
debug <- "I'm here"

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "output.png"
  args[3] = 0
}

data<-read.delim(args[1],header=T,sep="\t")
debug
if (args[3]==0){
#Grafico Normal -> opt=0
p<-ggplot(data, aes(x = reorder(Feature, Counts), y = Counts, fill = Feature)) +
 geom_bar(stat = "identity", show.legend = FALSE) +  # Eliminar leyenda
 coord_flip() +  # Voltear el gráfico para mejor lectura
 theme_bw(base_size = 14) +geom_text(aes(label = Counts), hjust = -0.1, size = 3)+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
 panel.background = element_blank()) +ggtitle("")+ xlab("") +ylab("")+scale_y_continuous(limits=c(0,170000),breaks = seq(10000,170000,20000))
}else{
#Grafico Volteado -> opt=1
p<-ggplot(data, aes(x = reorder(Feature, Counts), y = Counts, fill = Feature)) +
 geom_bar(stat = "identity", show.legend = FALSE) +  # Eliminar leyenda
 coord_flip() +  # Voltear el gráfico para mejor lectura
 theme_bw(base_size = 14) +geom_text(aes(label = Counts), hjust = 1.1, size = 3)+
 theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
 panel.background = element_blank()) +ggtitle("")+ xlab("") +ylab("")+scale_y_reverse(limits=c(170000,0),breaks = seq(10000,170000,20000))+scale_x_discrete(position = "top")
}

ggsave(p,filename=args[2],dpi=1000,width=12,height=8)
