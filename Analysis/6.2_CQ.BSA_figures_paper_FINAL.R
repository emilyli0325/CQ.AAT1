
library(magrittr)
library(dplyr)
library("QTLseqr")
library("ggpubr")

library(reshape2)
library(ggplot2)
library(ggridges)
library(doBy)

format_genomic <- function(...) {
      function(x) {
            limits <- c(1e0,   1e3, 1e6)
            #prefix <- c("","Kb","Mb")
            # Vector with array indices according to position in intervals
            i <- findInterval(abs(x), limits)
            # Set prefix to " " for very small values < 1e-24
            i <- ifelse(i==0, which(limits == 1e0), i)
            paste(format(round(x/limits[i], 1),
                         trim=TRUE, scientific=FALSE, ...)
                #  ,prefix[i]
            )
      }
}

################################################################
# NF54xNHP4026
################################################################

setwd("D:/Dropbox (TX Biomed)/Emily/4.P01/5.5.BSA5_NF54xNHP4026_all/QTL/Paper.plot")
setwd("C:/Users/xli.TXBIOMED/Dropbox (TX Biomed)/Emily/4.P01/5.5.BSA5_NF54xNHP4026_all/QTL/Paper.plot")


#### main figures

# allele frequency plots
refFre.AD.BSA5.LC <- read.csv(file = "C:/Users/xli.TXBIOMED/Dropbox (TX Biomed)/Emily/4.P01/5.5.BSA5_NF54xNHP4026_all/QTL/BSA5.ALL.AF.filter.tricube_coregenome.csv",header = TRUE)
refFre.AD.BSA5.LC$CHROM <- gsub("\\_v3|Pf3D7\\_0|Pf3D7\\_","",refFre.AD.BSA5.LC$CHROM)
refFre.AD.BSA5.LC$CHROM <- as.numeric(as.character(refFre.AD.BSA5.LC$CHROM))
p0 <- ggplot(data=refFre.AD.BSA5.LC)+scale_x_continuous(breaks=seq(from=0,to=max(refFre.AD.BSA5.LC$POS),by=10^(floor(log10(max(refFre.AD.BSA5.LC$POS))))),labels=format_genomic())+ylim(0,1)+facet_grid(~CHROM,scales="free_x",space="free_x") +theme_classic()


CQ.M13.48h.day4 <- p0 + ylab("3D7 allele frequency")+
     geom_line(aes_string(x = "POS", y = "AB_BC_0038.tricube"),color = "#000000",size=.6) + # before, black
     geom_line(aes_string(x = "POS", y = "AB_BC_0045.tricube"),color = "#009E73",size=.3) +  # d4.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0046.tricube"),color = "#009E73",size=.3) +  # d4.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0051.tricube"),color = "#56B4E9",size=.3) + #d4.50nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0052.tricube"),color = "#56B4E9",size=.3) + #d4.50nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0049.tricube"),color = "#E69F00",size=.3) + #d4.100nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0050.tricube"),color = "#E69F00",size=.3) + #d4.100nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0047.tricube"),color = "#D55E00",size=.3) + #d4.250nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0048.tricube"),color = "#D55E00",size=.3) #d4.250nM
	 
CQ.M23.48h.day4 <- p0 + ylab("3D7 allele frequency")+
     geom_line(aes_string(x = "POS", y = "AB_BC_0040.tricube"),color = "#000000",size=.6) + # before
     geom_line(aes_string(x = "POS", y = "AB_BC_0085.tricube"),color = "#009E73",size=.3) +  # d4.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0086.tricube"),color = "#009E73",size=.3) +  # d4.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0091.tricube"),color = "#56B4E9",size=.3) + #d4.50nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0092.tricube"),color = "#56B4E9",size=.3) + #d4.50nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0089.tricube"),color = "#E69F00",size=.3) + #d4.100nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0090.tricube"),color = "#E69F00",size=.3) + #d4.100nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0087.tricube"),color = "#D55E00",size=.3) + #d4.250nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0088.tricube"),color = "#D55E00",size=.3) #d4.250nM

CQ.M13.96h.day5 <- p0 + ylab("3D7 allele frequency")+
     geom_line(aes_string(x = "POS", y = "AB_BC_0038.tricube"),color = "#000000",size=.6) + # before
     geom_line(aes_string(x = "POS", y = "AB_BC_0165.tricube"),color = "#009E73",size=.3) +  # d5.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0166.tricube"),color = "#009E73",size=.3) +  # d5.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0171.tricube"),color = "#56B4E9",size=.3) + #d5.50nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0172.tricube"),color = "#56B4E9",size=.3) + #d5.50nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0169.tricube"),color = "#E69F00",size=.3) + #d5.100nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0170.tricube"),color = "#E69F00",size=.3) + #d5.100nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0167.tricube"),color = "#D55E00",size=.3) + #d5.250nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0168.tricube"),color = "#D55E00",size=.3) #d5.250nM

CQ.M23.96h.day5 <- p0 + ylab("3D7 allele frequency")+
     geom_line(aes_string(x = "POS", y = "AB_BC_0040.tricube"),color = "#000000",size=.6) + # before
     geom_line(aes_string(x = "POS", y = "AB_BC_0205.tricube"),color = "#009E73",size=.3) +  # d5.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0206.tricube"),color = "#009E73",size=.3) +  # d5.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0211.tricube"),color = "#56B4E9",size=.3) + #d5.50nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0212.tricube"),color = "#56B4E9",size=.3) + #d5.50nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0209.tricube"),color = "#E69F00",size=.3) + #d5.100nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0210.tricube"),color = "#E69F00",size=.3) + #d5.100nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0207.tricube"),color = "#D55E00",size=.3) + #d5.250nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0208.tricube"),color = "#D55E00",size=.3)  #d5.250nM	 


	 
CQ.M13.250nM <- p0 + ylab("3D7 allele frequency")+	 
     geom_line(aes_string(x = "POS", y = "AB_BC_0038.tricube"),color = "#000000",size=.6) + # before, black	
     geom_line(aes_string(x = "POS", y = "AB_BC_0045.tricube"),color = "#009E73",size=.3) +  # d4.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0046.tricube"),color = "#009E73",size=.3) +  # d4.control	 
	 geom_line(aes_string(x = "POS", y = "AB_BC_0047.tricube"),color = "#D55E00",size=.3) + #d4.250nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0048.tricube"),color = "#D55E00",size=.3) + #d4.250nM
     geom_line(aes_string(x = "POS", y = "AB_BC_0165.tricube"),color = "#56B4E9",size=.3) +  # d5.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0166.tricube"),color = "#56B4E9",size=.3) +  # d5.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0167.tricube"),color = "#0072B2",size=.3) + #d5.250nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0168.tricube"),color = "#0072B2",size=.3) #d5.250nM	 
	 
CQ.M23.250nM <- p0 + ylab("3D7 allele frequency")+
     geom_line(aes_string(x = "POS", y = "AB_BC_0040.tricube"),color = "#000000",size=.6) + # before
     geom_line(aes_string(x = "POS", y = "AB_BC_0085.tricube"),color = "#009E73",size=.3) +  # d4.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0086.tricube"),color = "#009E73",size=.3) +  # d4.control	 
	 geom_line(aes_string(x = "POS", y = "AB_BC_0087.tricube"),color = "#D55E00",size=.3) + #d4.250nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0088.tricube"),color = "#D55E00",size=.3) + #d4.250nM	 
     geom_line(aes_string(x = "POS", y = "AB_BC_0205.tricube"),color = "#56B4E9",size=.3) +  # d5.control
	 geom_line(aes_string(x = "POS", y = "AB_BC_0206.tricube"),color = "#56B4E9",size=.3) +  # d5.control	 
	 geom_line(aes_string(x = "POS", y = "AB_BC_0207.tricube"),color = "#0072B2",size=.3) + #d5.250nM
	 geom_line(aes_string(x = "POS", y = "AB_BC_0208.tricube"),color = "#0072B2",size=.3)  #d5.250nM		 
	 
	 
	 
# G plot 48h
QTL <- read.csv(file = "C:/Users/xli.TXBIOMED/Dropbox (TX Biomed)/Emily/4.P01/5.5.BSA5_NF54xNHP4026_all/QTL/QTLs.CQ/G_CQ.csv",header = TRUE)
QTL$CHROM <- gsub("\\_v3|Pf3D7\\_0|Pf3D7\\_","",QTL$CHROM)
QTL$CHROM <- as.numeric(as.character(QTL$CHROM))
QTL.day4 <- QTL[,c("M13.48h.d4.50.rep1","M13.48h.d4.50.rep2","M13.48h.d4.100.rep1","M13.48h.d4.100.rep2","M13.48h.d4.250.rep1","M13.48h.d4.250.rep2","M23.48h.d4.50.rep1","M23.48h.d4.50.rep2","M23.48h.d4.100.rep1","M23.48h.d4.100.rep2","M23.48h.d4.250.rep1","M23.48h.d4.250.rep2")]
QTL.day4 <- cbind(QTL[,1:5], QTL.day4)
QTL.filt <- QTL.day4[complete.cases(QTL.day4), ]

p0.QTL <- ggplot(data=QTL.filt)+scale_x_continuous(breaks=seq(from=0,to=max(QTL.filt$POS),by=10^(floor(log10(max(QTL.filt$POS))))),labels=format_genomic())+facet_grid(~CHROM,scales="free_x",space="free_x") +theme_classic()

p.M13.48h.d4 <- p0.QTL + ylab("G prime")+ ylim(0, 150)+
     geom_hline(yintercept = 20,color = "black",size = 0.6, linetype = "dashed") +
     geom_line(aes_string(x = "POS", y = "M13.48h.d4.50.rep1"),color = "#56B4E9",size=0.3) +
     geom_line(aes_string(x = "POS", y = "M13.48h.d4.50.rep2"),color = "#56B4E9",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.100.rep1"),color = "#E69F00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.100.rep2"),color = "#E69F00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.250.rep1"),color = "#D55E00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.250.rep2"),color = "#D55E00",size=0.3) 

p.M23.48h.d4 <- p0.QTL + ylab("G prime")+ ylim(0, 150)+
     geom_hline(yintercept = 20,color = "black",size = 0.6, linetype = "dashed") +
     geom_line(aes_string(x = "POS", y = "M23.48h.d4.50.rep1"),color = "#56B4E9",size=0.3) +
     geom_line(aes_string(x = "POS", y = "M23.48h.d4.50.rep2"),color = "#56B4E9",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.100.rep1"),color = "#E69F00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.100.rep2"),color = "#E69F00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.250.rep1"),color = "#D55E00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.250.rep2"),color = "#D55E00",size=0.3) 


	 
	 
	 
# G plot 96h	  
QTL.day5 <- QTL[,c("M13.96h.d5.50.rep1","M13.96h.d5.50.rep2","M13.96h.d5.100.rep1","M13.96h.d5.100.rep2","M13.96h.d5.250.rep1","M13.96h.d5.250.rep2","M23.96h.d5.50.rep1","M23.96h.d5.50.rep2","M23.96h.d5.100.rep1","M23.96h.d5.100.rep2","M23.96h.d5.250.rep1","M23.96h.d5.250.rep2")]
QTL.day5 <- cbind(QTL[,1:5], QTL.day5)
QTL.filt.d5 <- QTL.day5[complete.cases(QTL.day5), ]

p0.QTL.d5 <- ggplot(data=QTL.filt.d5)+scale_x_continuous(breaks=seq(from=0,to=max(QTL.filt.d5$POS),by=10^(floor(log10(max(QTL.filt.d5$POS))))),labels=format_genomic())+facet_grid(~CHROM,scales="free_x",space="free_x") +theme_classic()
	 
p.M13.96h.d5 <- p0.QTL.d5 + ylab("G prime")+ ylim(0, 150)+
     geom_hline(yintercept = 20,color = "black",size = 0.6, linetype = "dashed") +
     geom_line(aes_string(x = "POS", y = "M13.96h.d5.50.rep1"),color = "#56B4E9",size=0.3) +
     geom_line(aes_string(x = "POS", y = "M13.96h.d5.50.rep2"),color = "#56B4E9",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.100.rep1"),color = "#E69F00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.100.rep2"),color = "#E69F00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.250.rep1"),color = "#D55E00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.250.rep2"),color = "#D55E00",size=0.3) 

p.M23.96h.d5 <- p0.QTL.d5 + ylab("G prime")+ ylim(0, 150)+
     geom_hline(yintercept = 20,color = "black",size = 0.6, linetype = "dashed") +
     geom_line(aes_string(x = "POS", y = "M23.96h.d5.50.rep1"),color = "#56B4E9",size=0.3) +
     geom_line(aes_string(x = "POS", y = "M23.96h.d5.50.rep2"),color = "#56B4E9",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.100.rep1"),color = "#E69F00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.100.rep2"),color = "#E69F00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.250.rep1"),color = "#D55E00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.250.rep2"),color = "#D55E00",size=0.3) 	 
	 
	 
	 
pdf('CQ.BSA.NF54xNHP4026.48h.day4.pdf', width=6.5, height=10)
ggarrange(CQ.M13.48h.day4, CQ.M23.48h.day4, p.M13.48h.d4, p.M23.48h.d4,
          labels = c("","","",""),
          ncol = 1, nrow = 4)
dev.off()


pdf('CQ.BSA.NF54xNHP4026.96h.day5.pdf', width=6.5, height=10)
ggarrange(CQ.M13.96h.day5, CQ.M23.96h.day5, p.M13.96h.d5, p.M23.96h.d5,
          labels = c("","","",""),
          ncol = 1, nrow = 4)
dev.off()



#### 250nM only 
  
QTL.250nM <- QTL[,c("M13.48h.d4.250.rep1","M13.48h.d4.250.rep2","M23.48h.d4.250.rep1","M23.48h.d4.250.rep2","M13.96h.d5.250.rep1","M13.96h.d5.250.rep2","M23.96h.d5.250.rep1","M23.96h.d5.250.rep2")]
QTL.250nM <- cbind(QTL[,1:5], QTL.250nM)
QTL.filt.250nM <- QTL.250nM[complete.cases(QTL.250nM), ]

p0.QTL.250nM <- ggplot(data=QTL.filt.250nM)+scale_x_continuous(breaks=seq(from=0,to=max(QTL.filt.250nM$POS),by=10^(floor(log10(max(QTL.filt.250nM$POS))))),labels=format_genomic())+facet_grid(~CHROM,scales="free_x",space="free_x") +theme_classic()
	 
p.M13.250nM <- p0.QTL.250nM + ylab("G prime")+ ylim(0, 150)+
     geom_hline(yintercept = 20,color = "black",size = 0.6, linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.250.rep1"),color = "#D55E00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.250.rep2"),color = "#D55E00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.250.rep1"),color = "#0072B2",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.250.rep2"),color = "#0072B2",size=0.3) 

p.M23.250nM <- p0.QTL.250nM + ylab("G prime")+ ylim(0, 150)+
     geom_hline(yintercept = 20,color = "black",size = 0.6, linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.250.rep1"),color = "#D55E00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.250.rep2"),color = "#D55E00",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.250.rep1"),color = "#0072B2",size=0.3) +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.250.rep2"),color = "#0072B2",size=0.3) 	 
	 


pdf('CQ.BSA.NF54xNHP4026.250nM.pdf', width=6.5, height=10)
ggarrange(CQ.M13.250nM, CQ.M23.250nM, p.M13.250nM, p.M23.250nM,
          labels = c("","","",""),
          ncol = 1, nrow = 4)
dev.off()


################################################################
# fine mapping of chr6 QTL

QTL.chr6 <- read.csv(file = "C:/Users/xli.TXBIOMED/Dropbox (TX Biomed)/Emily/4.P01/5.5.BSA5_NF54xNHP4026_all/QTL/QTLs.CQ/G_CQ_chr6.csv",header = TRUE)

QTL.chr6 <- read.csv(file = "D:/Dropbox (TX Biomed)/Emily/4.P01/5.5.BSA5_NF54xNHP4026_all/QTL/QTLs.CQ/G_CQ_chr6.csv",header = TRUE)

QTL.chr6$CHROM <- gsub("\\_v3|Pf3D7\\_0|Pf3D7\\_","",QTL.chr6$CHROM)
QTL.chr6$CHROM <- as.numeric(as.character(QTL.chr6$CHROM))

QTL.chr6.day4 <- QTL.chr6[,c("M13.48h.d4.250.rep1","M13.48h.d4.250.rep2","M13.96h.d5.250.rep1","M13.96h.d5.250.rep2","M23.48h.d4.250.rep1","M23.48h.d4.250.rep2","M23.96h.d5.250.rep1","M23.96h.d5.250.rep2")]
QTL.chr6.day4 <- cbind(QTL.chr6[,1:5], QTL.chr6.day4)
QTL.chr6.filt <- QTL.chr6.day4[complete.cases(QTL.chr6.day4), ]


## chr6
p0.QTL.chr6.filt <- ggplot(data=QTL.chr6.filt)+
      scale_x_continuous(breaks=seq(from=0,to=max(QTL.chr6.filt$POS),by=10^(floor(log10(max(QTL.chr6.filt$POS))))),labels=format_genomic())+
	  facet_grid(~CHROM,scales="free_x",space="free_x") +
	  theme(plot.margin = margin(b = 10,l = 20,r = 20,unit = "pt")) + theme_bw()
	  
p.chr6 <- p0.QTL.chr6.filt + ylab(" G prime ")+
     geom_rect(data = QTL.chr6.filt, aes(xmin = 1012526, xmax = 1283256, ymin = -Inf, ymax = Inf), fill = "lightcyan2", alpha = 0.4)+
	 geom_hline(yintercept = 20,color = "black",size = 0.5, linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.250.rep1"),color = "#D55E00",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.250.rep2"),color = "#D55E00",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.250.rep1"),color = "#0072B2",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.250.rep2"),color = "#0072B2",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.250.rep1"),color = "#D55E00",size=0.6,linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.250.rep2"),color = "#D55E00",size=0.6,linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.250.rep1"),color = "#0072B2",size=0.6,linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.250.rep2"),color = "#0072B2",size=0.6,linetype = "dashed")

p.chr6

pdf('chr6.QTL.region.G.pdf', width=3, height=4)
p.chr6
dev.off()




p.chr6 <- p0.QTL.chr6.filt + ylab(" G prime ")+
     geom_rect(data = QTL.chr6.filt, aes(xmin = 1012526, xmax = 1283256, ymin = -Inf, ymax = Inf), fill = "lightcyan2", alpha = 0.4)+
	 geom_hline(yintercept = 20,color = "black",size = 0.5, linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.250.rep1"),color = "#0072B2",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.250.rep2"),color = "#0072B2",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.250.rep1"),color = "#D55E00",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.250.rep2"),color = "#D55E00",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.250.rep1"),color = "#0072B2",size=0.6,linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.250.rep2"),color = "#0072B2",size=0.6,linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.250.rep1"),color = "#D55E00",size=0.6,linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.250.rep2"),color = "#D55E00",size=0.6,linetype = "dashed")

p.chr6

pdf('chr6.QTL.region.G_version2.pdf', width=3, height=4)
p.chr6
dev.off()




p.chr6 <- p0.QTL.chr6.filt + ylab(" G prime ")+
     geom_rect(data = QTL.chr6.filt, aes(xmin = 1012526, xmax = 1283256, ymin = -Inf, ymax = Inf), fill = "lightcyan2", alpha = 0.4)+
	 geom_hline(yintercept = 20,color = "black",size = 0.5, linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.250.rep1"),color = "#D55E00",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M13.48h.d4.250.rep2"),color = "#D55E00",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.250.rep1"),color = "#0072B2",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M13.96h.d5.250.rep2"),color = "#0072B2",size=0.6) +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.250.rep1"),color = "#D55E00",size=0.6,linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M23.48h.d4.250.rep2"),color = "#D55E00",size=0.6,linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.250.rep1"),color = "#0072B2",size=0.6,linetype = "dashed") +
	 geom_line(aes_string(x = "POS", y = "M23.96h.d5.250.rep2"),color = "#0072B2",size=0.6,linetype = "dashed") +
	 geom_vline(xintercept = 1213102, linetype="dashed", color = "red", size=1)  

p.chr6

pdf('chr6.QTL.region.G_version3.pdf', width=3, height=4)
p.chr6
dev.off()













