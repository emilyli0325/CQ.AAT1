######## BSA CQ P, G
setwd("D:/Dropbox (TX Biomed)/Emily/P01/5.5.BSA5_NF54xNHP4026_all/QTL")
library(magrittr)
library(dplyr)
library("QTLseqr")
library("ggpubr")
library(ggpubr)
################################################################################################################
################################################################################################################

# M13 48h day4 - control
d1.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0045",lowBulk = "AB_BC_0051")
d1.1.filt <-filterSNPs(SNPset = d1.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d1.1.filt <- runQTLseqAnalysis(SNPset = d1.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d1.1.filt <- runGprimeAnalysis(d1.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d1.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d1.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day4.50nM-control.rep1.csv",row.names=FALSE)

d1.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0046",lowBulk = "AB_BC_0052")
d1.2.filt <-filterSNPs(SNPset = d1.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d1.2.filt <- runQTLseqAnalysis(SNPset = d1.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d1.2.filt <- runGprimeAnalysis(d1.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d1.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d1.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day4.50nM-control.rep2.csv",row.names=FALSE)


d1.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0045",lowBulk = "AB_BC_0049")
d1.3.filt <-filterSNPs(SNPset = d1.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d1.3.filt <- runQTLseqAnalysis(SNPset = d1.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d1.3.filt <- runGprimeAnalysis(d1.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d1.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d1.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day4.100nM-control.rep1.csv",row.names=FALSE)


d1.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0046",lowBulk = "AB_BC_0050")
d1.4.filt <-filterSNPs(SNPset = d1.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d1.4.filt <- runQTLseqAnalysis(SNPset = d1.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d1.4.filt <- runGprimeAnalysis(d1.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d1.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d1.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day4.100nM-control.rep2.csv",row.names=FALSE)

d1.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0045",lowBulk = "AB_BC_0047")
d1.5.filt <-filterSNPs(SNPset = d1.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d1.5.filt <- runQTLseqAnalysis(SNPset = d1.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d1.5.filt <- runGprimeAnalysis(d1.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d1.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d1.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day4.250nM-control.rep1.csv",row.names=FALSE)

d1.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0046",lowBulk = "AB_BC_0048")
d1.6.filt <-filterSNPs(SNPset = d1.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d1.6.filt <- runQTLseqAnalysis(SNPset = d1.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d1.6.filt <- runGprimeAnalysis(d1.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d1.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d1.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day4.250nM-control.rep2.csv",row.names=FALSE)



################################################################################################################
################################################################################################################

# M13 48h day7 - control
d2.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0349",lowBulk = "AB_BC_0355")
d2.1.filt <-filterSNPs(SNPset = d2.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d2.1.filt <- runQTLseqAnalysis(SNPset = d2.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d2.1.filt <- runGprimeAnalysis(d2.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d2.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d2.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day7.50nM-control.rep1.csv",row.names=FALSE)
d2.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0350",lowBulk = "AB_BC_0356")
d2.2.filt <-filterSNPs(SNPset = d2.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d2.2.filt <- runQTLseqAnalysis(SNPset = d2.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d2.2.filt <- runGprimeAnalysis(d2.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d2.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d2.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day7.50nM-control.rep2.csv",row.names=FALSE)


d2.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0349",lowBulk = "AB_BC_0353")
d2.3.filt <-filterSNPs(SNPset = d2.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d2.3.filt <- runQTLseqAnalysis(SNPset = d2.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d2.3.filt <- runGprimeAnalysis(d2.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d2.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d2.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day7.100nM-control.rep1.csv",row.names=FALSE)
d2.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0350",lowBulk = "AB_BC_0354")
d2.4.filt <-filterSNPs(SNPset = d2.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d2.4.filt <- runQTLseqAnalysis(SNPset = d2.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d2.4.filt <- runGprimeAnalysis(d2.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d2.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d2.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day7.100nM-control.rep2.csv",row.names=FALSE)


d2.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0349",lowBulk = "AB_BC_0351")
d2.5.filt <-filterSNPs(SNPset = d2.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d2.5.filt <- runQTLseqAnalysis(SNPset = d2.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d2.5.filt <- runGprimeAnalysis(d2.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d2.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d2.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day7.250nM-control.rep1.csv",row.names=FALSE)
d2.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0350",lowBulk = "AB_BC_0352")
d2.6.filt <-filterSNPs(SNPset = d2.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d2.6.filt <- runQTLseqAnalysis(SNPset = d2.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d2.6.filt <- runGprimeAnalysis(d2.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d2.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d2.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.48h.day7.250nM-control.rep2.csv",row.names=FALSE)

################################################################################################################
################################################################################################################

# M13 96h day5 - control
d3.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0165",lowBulk = "AB_BC_0171")
d3.1.filt <-filterSNPs(SNPset = d3.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d3.1.filt <- runQTLseqAnalysis(SNPset = d3.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d3.1.filt <- runGprimeAnalysis(d3.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d3.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d3.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day5.50nM-control.rep1.csv",row.names=FALSE)
d3.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0166",lowBulk = "AB_BC_0172")
d3.2.filt <-filterSNPs(SNPset = d3.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d3.2.filt <- runQTLseqAnalysis(SNPset = d3.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d3.2.filt <- runGprimeAnalysis(d3.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d3.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d3.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day5.50nM-control.rep2.csv",row.names=FALSE)


d3.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0165",lowBulk = "AB_BC_0169")
d3.3.filt <-filterSNPs(SNPset = d3.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d3.3.filt <- runQTLseqAnalysis(SNPset = d3.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d3.3.filt <- runGprimeAnalysis(d3.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d3.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d3.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day5.100nM-control.rep1.csv",row.names=FALSE)
d3.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0166",lowBulk = "AB_BC_0170")
d3.4.filt <-filterSNPs(SNPset = d3.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d3.4.filt <- runQTLseqAnalysis(SNPset = d3.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d3.4.filt <- runGprimeAnalysis(d3.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d3.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d3.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day5.100nM-control.rep2.csv",row.names=FALSE)

d3.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0165",lowBulk = "AB_BC_0167")
d3.5.filt <-filterSNPs(SNPset = d3.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d3.5.filt <- runQTLseqAnalysis(SNPset = d3.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d3.5.filt <- runGprimeAnalysis(d3.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d3.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d3.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day5.250nM-control.rep1.csv",row.names=FALSE)
d3.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0166",lowBulk = "AB_BC_0168")
d3.6.filt <-filterSNPs(SNPset = d3.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d3.6.filt <- runQTLseqAnalysis(SNPset = d3.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d3.6.filt <- runGprimeAnalysis(d3.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d3.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d3.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day5.250nM-control.rep2.csv",row.names=FALSE)

################################################################################################################
################################################################################################################

# M13 96h day10 - control
d4.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0469",lowBulk = "AB_BC_0475")
d4.1.filt <-filterSNPs(SNPset = d4.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d4.1.filt <- runQTLseqAnalysis(SNPset = d4.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d4.1.filt <- runGprimeAnalysis(d4.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d4.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d4.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day10.50nM-control.rep1.csv",row.names=FALSE)
d4.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0470",lowBulk = "AB_BC_0476")
d4.2.filt <-filterSNPs(SNPset = d4.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d4.2.filt <- runQTLseqAnalysis(SNPset = d4.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d4.2.filt <- runGprimeAnalysis(d4.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d4.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d4.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day10.50nM-control.rep2.csv",row.names=FALSE)


d4.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0469",lowBulk = "AB_BC_0473")
d4.3.filt <-filterSNPs(SNPset = d4.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d4.3.filt <- runQTLseqAnalysis(SNPset = d4.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d4.3.filt <- runGprimeAnalysis(d4.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d4.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d4.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day10.100nM-control.rep1.csv",row.names=FALSE)
d4.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0470",lowBulk = "AB_BC_0474")
d4.4.filt <-filterSNPs(SNPset = d4.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d4.4.filt <- runQTLseqAnalysis(SNPset = d4.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d4.4.filt <- runGprimeAnalysis(d4.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d4.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d4.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day10.100nM-control.rep2.csv",row.names=FALSE)

d4.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0469",lowBulk = "AB_BC_0471")
d4.5.filt <-filterSNPs(SNPset = d4.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d4.5.filt <- runQTLseqAnalysis(SNPset = d4.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d4.5.filt <- runGprimeAnalysis(d4.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d4.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d4.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day10.250nM-control.rep1.csv",row.names=FALSE)
d4.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0470",lowBulk = "AB_BC_0472")
d4.6.filt <-filterSNPs(SNPset = d4.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d4.6.filt <- runQTLseqAnalysis(SNPset = d4.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d4.6.filt <- runGprimeAnalysis(d4.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d4.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d4.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M13.96h.day10.250nM-control.rep2.csv",row.names=FALSE)


################################################################################################################
################################################################################################################

# M16 48h day4 - control
d5.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0065",lowBulk = "AB_BC_0071")
d5.1.filt <-filterSNPs(SNPset = d5.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d5.1.filt <- runQTLseqAnalysis(SNPset = d5.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d5.1.filt <- runGprimeAnalysis(d5.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d5.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d5.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day4.50nM-control.rep1.csv",row.names=FALSE)
d5.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0066",lowBulk = "AB_BC_0072")
d5.2.filt <-filterSNPs(SNPset = d5.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d5.2.filt <- runQTLseqAnalysis(SNPset = d5.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d5.2.filt <- runGprimeAnalysis(d5.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d5.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d5.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day4.50nM-control.rep2.csv",row.names=FALSE)

d5.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0065",lowBulk = "AB_BC_0069")
d5.3.filt <-filterSNPs(SNPset = d5.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d5.3.filt <- runQTLseqAnalysis(SNPset = d5.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d5.3.filt <- runGprimeAnalysis(d5.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d5.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d5.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day4.100nM-control.rep1.csv",row.names=FALSE)
d5.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0066",lowBulk = "AB_BC_0070")
d5.4.filt <-filterSNPs(SNPset = d5.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d5.4.filt <- runQTLseqAnalysis(SNPset = d5.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d5.4.filt <- runGprimeAnalysis(d5.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d5.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d5.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day4.100nM-control.rep2.csv",row.names=FALSE)

d5.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0065",lowBulk = "AB_BC_0067")
d5.5.filt <-filterSNPs(SNPset = d5.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d5.5.filt <- runQTLseqAnalysis(SNPset = d5.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d5.5.filt <- runGprimeAnalysis(d5.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d5.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d5.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day4.250nM-control.rep1.csv",row.names=FALSE)
d5.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0066",lowBulk = "AB_BC_0068")
d5.6.filt <-filterSNPs(SNPset = d5.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d5.6.filt <- runQTLseqAnalysis(SNPset = d5.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d5.6.filt <- runGprimeAnalysis(d5.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d5.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d5.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day4.250nM-control.rep2.csv",row.names=FALSE)

################################################################################################################
################################################################################################################

# M16 48h day7 - control
d6.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0369",lowBulk = "AB_BC_0375")
d6.1.filt <-filterSNPs(SNPset = d6.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d6.1.filt <- runQTLseqAnalysis(SNPset = d6.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d6.1.filt <- runGprimeAnalysis(d6.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d6.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d6.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day7.50nM-control.rep1.csv",row.names=FALSE)
d6.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0370",lowBulk = "AB_BC_0376")
d6.2.filt <-filterSNPs(SNPset = d6.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d6.2.filt <- runQTLseqAnalysis(SNPset = d6.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d6.2.filt <- runGprimeAnalysis(d6.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d6.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d6.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day7.50nM-control.rep2.csv",row.names=FALSE)


d6.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0369",lowBulk = "AB_BC_0373")
d6.3.filt <-filterSNPs(SNPset = d6.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d6.3.filt <- runQTLseqAnalysis(SNPset = d6.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d6.3.filt <- runGprimeAnalysis(d6.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d6.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d6.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day7.100nM-control.rep1.csv",row.names=FALSE)
d6.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0370",lowBulk = "AB_BC_0374")
d6.4.filt <-filterSNPs(SNPset = d6.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d6.4.filt <- runQTLseqAnalysis(SNPset = d6.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d6.4.filt <- runGprimeAnalysis(d6.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d6.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d6.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day7.100nM-control.rep2.csv",row.names=FALSE)

d6.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0369",lowBulk = "AB_BC_0371")
d6.5.filt <-filterSNPs(SNPset = d6.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d6.5.filt <- runQTLseqAnalysis(SNPset = d6.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d6.5.filt <- runGprimeAnalysis(d6.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d6.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d6.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day7.250nM-control.rep1.csv",row.names=FALSE)
d6.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0370",lowBulk = "AB_BC_0372")
d6.6.filt <-filterSNPs(SNPset = d6.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d6.6.filt <- runQTLseqAnalysis(SNPset = d6.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d6.6.filt <- runGprimeAnalysis(d6.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d6.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d6.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.48h.day7.250nM-control.rep2.csv",row.names=FALSE)


################################################################################################################
################################################################################################################

# M16 96h day5 - control
d7.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0185",lowBulk = "AB_BC_0191")
d7.1.filt <-filterSNPs(SNPset = d7.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d7.1.filt <- runQTLseqAnalysis(SNPset = d7.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d7.1.filt <- runGprimeAnalysis(d7.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d7.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d7.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day5.50nM-control.rep1.csv",row.names=FALSE)
d7.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0186",lowBulk = "AB_BC_0192")
d7.2.filt <-filterSNPs(SNPset = d7.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d7.2.filt <- runQTLseqAnalysis(SNPset = d7.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d7.2.filt <- runGprimeAnalysis(d7.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d7.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d7.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day5.50nM-control.rep2.csv",row.names=FALSE)


d7.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0185",lowBulk = "AB_BC_0189")
d7.3.filt <-filterSNPs(SNPset = d7.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d7.3.filt <- runQTLseqAnalysis(SNPset = d7.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d7.3.filt <- runGprimeAnalysis(d7.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d7.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d7.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day5.100nM-control.rep1.csv",row.names=FALSE)
d7.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0186",lowBulk = "AB_BC_0190")
d7.4.filt <-filterSNPs(SNPset = d7.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d7.4.filt <- runQTLseqAnalysis(SNPset = d7.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d7.4.filt <- runGprimeAnalysis(d7.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d7.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d7.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day5.100nM-control.rep2.csv",row.names=FALSE)

d7.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0185",lowBulk = "AB_BC_0187")
d7.5.filt <-filterSNPs(SNPset = d7.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d7.5.filt <- runQTLseqAnalysis(SNPset = d7.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d7.5.filt <- runGprimeAnalysis(d7.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d7.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d7.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day5.250nM-control.rep1.csv",row.names=FALSE)
d7.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0186",lowBulk = "AB_BC_0188")
d7.6.filt <-filterSNPs(SNPset = d7.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d7.6.filt <- runQTLseqAnalysis(SNPset = d7.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d7.6.filt <- runGprimeAnalysis(d7.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d7.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d7.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day5.250nM-control.rep2.csv",row.names=FALSE)


################################################################################################################
################################################################################################################

# M16 96h day10 - control
d8.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0489",lowBulk = "AB_BC_0495")
d8.1.filt <-filterSNPs(SNPset = d8.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d8.1.filt <- runQTLseqAnalysis(SNPset = d8.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d8.1.filt <- runGprimeAnalysis(d8.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d8.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d8.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day10.50nM-control.rep1.csv",row.names=FALSE)
d8.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0490",lowBulk = "AB_BC_0496")
d8.2.filt <-filterSNPs(SNPset = d8.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d8.2.filt <- runQTLseqAnalysis(SNPset = d8.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d8.2.filt <- runGprimeAnalysis(d8.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d8.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d8.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day10.50nM-control.rep2.csv",row.names=FALSE)


d8.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0489",lowBulk = "AB_BC_0493")
d8.3.filt <-filterSNPs(SNPset = d8.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d8.3.filt <- runQTLseqAnalysis(SNPset = d8.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d8.3.filt <- runGprimeAnalysis(d8.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d8.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d8.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day10.100nM-control.rep1.csv",row.names=FALSE)
d8.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0490",lowBulk = "AB_BC_0494")
d8.4.filt <-filterSNPs(SNPset = d8.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d8.4.filt <- runQTLseqAnalysis(SNPset = d8.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d8.4.filt <- runGprimeAnalysis(d8.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d8.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d8.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day10.100nM-control.rep2.csv",row.names=FALSE)

d8.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0489",lowBulk = "AB_BC_0491")
d8.5.filt <-filterSNPs(SNPset = d8.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d8.5.filt <- runQTLseqAnalysis(SNPset = d8.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d8.5.filt <- runGprimeAnalysis(d8.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d8.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d8.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day10.250nM-control.rep1.csv",row.names=FALSE)
d8.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0490",lowBulk = "AB_BC_0492")
d8.6.filt <-filterSNPs(SNPset = d8.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d8.6.filt <- runQTLseqAnalysis(SNPset = d8.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d8.6.filt <- runGprimeAnalysis(d8.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d8.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d8.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M16.96h.day10.250nM-control.rep2.csv",row.names=FALSE)


################################################################################################################
################################################################################################################

# M23 48h day4 - control
d9.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0085",lowBulk = "AB_BC_0091")
d9.1.filt <-filterSNPs(SNPset = d9.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d9.1.filt <- runQTLseqAnalysis(SNPset = d9.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d9.1.filt <- runGprimeAnalysis(d9.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d9.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d9.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day4.50nM-control.rep1.csv",row.names=FALSE)
d9.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0086",lowBulk = "AB_BC_0092")
d9.2.filt <-filterSNPs(SNPset = d9.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d9.2.filt <- runQTLseqAnalysis(SNPset = d9.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d9.2.filt <- runGprimeAnalysis(d9.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d9.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d9.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day4.50nM-control.rep2.csv",row.names=FALSE)


d9.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0085",lowBulk = "AB_BC_0089")
d9.3.filt <-filterSNPs(SNPset = d9.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d9.3.filt <- runQTLseqAnalysis(SNPset = d9.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d9.3.filt <- runGprimeAnalysis(d9.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d9.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d9.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day4.100nM-control.rep1.csv",row.names=FALSE)
d9.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0086",lowBulk = "AB_BC_0090")
d9.4.filt <-filterSNPs(SNPset = d9.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d9.4.filt <- runQTLseqAnalysis(SNPset = d9.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d9.4.filt <- runGprimeAnalysis(d9.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d9.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d9.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day4.100nM-control.rep2.csv",row.names=FALSE)

d9.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0085",lowBulk = "AB_BC_0087")
d9.5.filt <-filterSNPs(SNPset = d9.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d9.5.filt <- runQTLseqAnalysis(SNPset = d9.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d9.5.filt <- runGprimeAnalysis(d9.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d9.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d9.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day4.250nM-control.rep1.csv",row.names=FALSE)
d9.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0086",lowBulk = "AB_BC_0088")
d9.6.filt <-filterSNPs(SNPset = d9.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d9.6.filt <- runQTLseqAnalysis(SNPset = d9.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d9.6.filt <- runGprimeAnalysis(d9.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d9.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d9.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day4.250nM-control.rep2.csv",row.names=FALSE)


################################################################################################################
################################################################################################################

# M23 48h day7 - control
d10.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0389",lowBulk = "AB_BC_0395")
d10.1.filt <-filterSNPs(SNPset = d10.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d10.1.filt <- runQTLseqAnalysis(SNPset = d10.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d10.1.filt <- runGprimeAnalysis(d10.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d10.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d10.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day7.50nM-control.rep1.csv",row.names=FALSE)
d10.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0390",lowBulk = "AB_BC_0396")
d10.2.filt <-filterSNPs(SNPset = d10.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d10.2.filt <- runQTLseqAnalysis(SNPset = d10.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d10.2.filt <- runGprimeAnalysis(d10.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d10.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d10.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day7.50nM-control.rep2.csv",row.names=FALSE)


d10.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0389",lowBulk = "AB_BC_0393")
d10.3.filt <-filterSNPs(SNPset = d10.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d10.3.filt <- runQTLseqAnalysis(SNPset = d10.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d10.3.filt <- runGprimeAnalysis(d10.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d10.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d10.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day7.100nM-control.rep1.csv",row.names=FALSE)
d10.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0390",lowBulk = "AB_BC_0394")
d10.4.filt <-filterSNPs(SNPset = d10.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d10.4.filt <- runQTLseqAnalysis(SNPset = d10.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d10.4.filt <- runGprimeAnalysis(d10.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d10.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d10.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day7.100nM-control.rep2.csv",row.names=FALSE)

d10.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0389",lowBulk = "AB_BC_0391")
d10.5.filt <-filterSNPs(SNPset = d10.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d10.5.filt <- runQTLseqAnalysis(SNPset = d10.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d10.5.filt <- runGprimeAnalysis(d10.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d10.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d10.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day7.250nM-control.rep1.csv",row.names=FALSE)
d10.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0390",lowBulk = "AB_BC_0392")
d10.6.filt <-filterSNPs(SNPset = d10.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d10.6.filt <- runQTLseqAnalysis(SNPset = d10.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d10.6.filt <- runGprimeAnalysis(d10.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d10.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d10.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.48h.day7.250nM-control.rep2.csv",row.names=FALSE)


################################################################################################################
################################################################################################################

# M23 96h day5 - control
d11.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0205",lowBulk = "AB_BC_0211")
d11.1.filt <-filterSNPs(SNPset = d11.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d11.1.filt <- runQTLseqAnalysis(SNPset = d11.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d11.1.filt <- runGprimeAnalysis(d11.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d11.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d11.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day5.50nM-control.rep1.csv",row.names=FALSE)
d11.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0206",lowBulk = "AB_BC_0212")
d11.2.filt <-filterSNPs(SNPset = d11.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d11.2.filt <- runQTLseqAnalysis(SNPset = d11.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d11.2.filt <- runGprimeAnalysis(d11.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d11.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d11.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day5.50nM-control.rep2.csv",row.names=FALSE)


d11.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0205",lowBulk = "AB_BC_0209")
d11.3.filt <-filterSNPs(SNPset = d11.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d11.3.filt <- runQTLseqAnalysis(SNPset = d11.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d11.3.filt <- runGprimeAnalysis(d11.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d11.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d11.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day5.100nM-control.rep1.csv",row.names=FALSE)
d11.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0206",lowBulk = "AB_BC_0210")
d11.4.filt <-filterSNPs(SNPset = d11.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d11.4.filt <- runQTLseqAnalysis(SNPset = d11.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d11.4.filt <- runGprimeAnalysis(d11.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d11.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d11.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day5.100nM-control.rep2.csv",row.names=FALSE)

d11.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0205",lowBulk = "AB_BC_0207")
d11.5.filt <-filterSNPs(SNPset = d11.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d11.5.filt <- runQTLseqAnalysis(SNPset = d11.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d11.5.filt <- runGprimeAnalysis(d11.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d11.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d11.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day5.250nM-control.rep1.csv",row.names=FALSE)
d11.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0206",lowBulk = "AB_BC_0208")
d11.6.filt <-filterSNPs(SNPset = d11.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d11.6.filt <- runQTLseqAnalysis(SNPset = d11.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d11.6.filt <- runGprimeAnalysis(d11.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d11.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d11.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day5.250nM-control.rep2.csv",row.names=FALSE)


################################################################################################################
################################################################################################################

# M23 96h day10 - control
d12.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0509",lowBulk = "AB_BC_0515")
d12.1.filt <-filterSNPs(SNPset = d12.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d12.1.filt <- runQTLseqAnalysis(SNPset = d12.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d12.1.filt <- runGprimeAnalysis(d12.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d12.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d12.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day10.50nM-control.rep1.csv",row.names=FALSE)
d12.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0510",lowBulk = "AB_BC_0516")
d12.2.filt <-filterSNPs(SNPset = d12.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d12.2.filt <- runQTLseqAnalysis(SNPset = d12.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d12.2.filt <- runGprimeAnalysis(d12.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d12.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d12.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day10.50nM-control.rep2.csv",row.names=FALSE)


d12.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0509",lowBulk = "AB_BC_0513")
d12.3.filt <-filterSNPs(SNPset = d12.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d12.3.filt <- runQTLseqAnalysis(SNPset = d12.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d12.3.filt <- runGprimeAnalysis(d12.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d12.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d12.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day10.100nM-control.rep1.csv",row.names=FALSE)
d12.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0510",lowBulk = "AB_BC_0514")
d12.4.filt <-filterSNPs(SNPset = d12.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d12.4.filt <- runQTLseqAnalysis(SNPset = d12.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d12.4.filt <- runGprimeAnalysis(d12.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d12.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d12.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day10.100nM-control.rep2.csv",row.names=FALSE)

d12.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0509",lowBulk = "AB_BC_0511")
d12.5.filt <-filterSNPs(SNPset = d12.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d12.5.filt <- runQTLseqAnalysis(SNPset = d12.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d12.5.filt <- runGprimeAnalysis(d12.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d12.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d12.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day10.250nM-control.rep1.csv",row.names=FALSE)
d12.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0510",lowBulk = "AB_BC_0512")
d12.6.filt <-filterSNPs(SNPset = d12.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d12.6.filt <- runQTLseqAnalysis(SNPset = d12.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d12.6.filt <- runGprimeAnalysis(d12.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d12.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d12.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M23.96h.day10.250nM-control.rep2.csv",row.names=FALSE)


################################################################################################################
################################################################################################################

# M26 48h day4 - control
d13.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0105",lowBulk = "AB_BC_0109")
d13.1.filt <-filterSNPs(SNPset = d13.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d13.1.filt <- runQTLseqAnalysis(SNPset = d13.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d13.1.filt <- runGprimeAnalysis(d13.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d13.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d13.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day4.50nM-control.rep1.csv",row.names=FALSE)
d13.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0106",lowBulk = "AB_BC_0110")
d13.2.filt <-filterSNPs(SNPset = d13.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d13.2.filt <- runQTLseqAnalysis(SNPset = d13.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d13.2.filt <- runGprimeAnalysis(d13.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d13.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d13.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day4.50nM-control.rep2.csv",row.names=FALSE)


d13.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0105",lowBulk = "AB_BC_0111")
d13.3.filt <-filterSNPs(SNPset = d13.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d13.3.filt <- runQTLseqAnalysis(SNPset = d13.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d13.3.filt <- runGprimeAnalysis(d13.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d13.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d13.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day4.100nM-control.rep1.csv",row.names=FALSE)
d13.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0106",lowBulk = "AB_BC_0112")
d13.4.filt <-filterSNPs(SNPset = d13.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d13.4.filt <- runQTLseqAnalysis(SNPset = d13.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d13.4.filt <- runGprimeAnalysis(d13.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d13.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d13.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day4.100nM-control.rep2.csv",row.names=FALSE)

d13.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0105",lowBulk = "AB_BC_0107")
d13.5.filt <-filterSNPs(SNPset = d13.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d13.5.filt <- runQTLseqAnalysis(SNPset = d13.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d13.5.filt <- runGprimeAnalysis(d13.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d13.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d13.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day4.250nM-control.rep1.csv",row.names=FALSE)
d13.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0106",lowBulk = "AB_BC_0108")
d13.6.filt <-filterSNPs(SNPset = d13.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d13.6.filt <- runQTLseqAnalysis(SNPset = d13.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d13.6.filt <- runGprimeAnalysis(d13.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d13.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d13.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day4.250nM-control.rep2.csv",row.names=FALSE)


################################################################################################################
################################################################################################################

# M26 48h day7 - control
d14.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0409",lowBulk = "AB_BC_0415")
d14.1.filt <-filterSNPs(SNPset = d14.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d14.1.filt <- runQTLseqAnalysis(SNPset = d14.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d14.1.filt <- runGprimeAnalysis(d14.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d14.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d14.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day7.50nM-control.rep1.csv",row.names=FALSE)
d14.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0410",lowBulk = "AB_BC_0416")
d14.2.filt <-filterSNPs(SNPset = d14.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d14.2.filt <- runQTLseqAnalysis(SNPset = d14.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d14.2.filt <- runGprimeAnalysis(d14.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d14.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d14.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day7.50nM-control.rep2.csv",row.names=FALSE)


d14.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0409",lowBulk = "AB_BC_0413")
d14.3.filt <-filterSNPs(SNPset = d14.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d14.3.filt <- runQTLseqAnalysis(SNPset = d14.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d14.3.filt <- runGprimeAnalysis(d14.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d14.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d14.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day7.100nM-control.rep1.csv",row.names=FALSE)
d14.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0410",lowBulk = "AB_BC_0414")
d14.4.filt <-filterSNPs(SNPset = d14.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d14.4.filt <- runQTLseqAnalysis(SNPset = d14.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d14.4.filt <- runGprimeAnalysis(d14.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d14.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d14.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day7.100nM-control.rep2.csv",row.names=FALSE)

d14.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0409",lowBulk = "AB_BC_0411")
d14.5.filt <-filterSNPs(SNPset = d14.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d14.5.filt <- runQTLseqAnalysis(SNPset = d14.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d14.5.filt <- runGprimeAnalysis(d14.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d14.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d14.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day7.250nM-control.rep1.csv",row.names=FALSE)
d14.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0410",lowBulk = "AB_BC_0412")
d14.6.filt <-filterSNPs(SNPset = d14.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d14.6.filt <- runQTLseqAnalysis(SNPset = d14.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d14.6.filt <- runGprimeAnalysis(d14.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d14.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d14.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.48h.day7.250nM-control.rep2.csv",row.names=FALSE)
################################################################################################################
################################################################################################################

# M26 96h day5 - control
d15.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0225",lowBulk = "AB_BC_0231")
d15.1.filt <-filterSNPs(SNPset = d15.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d15.1.filt <- runQTLseqAnalysis(SNPset = d15.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d15.1.filt <- runGprimeAnalysis(d15.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d15.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d15.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day5.50nM-control.rep1.csv",row.names=FALSE)
d15.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0226",lowBulk = "AB_BC_0232")
d15.2.filt <-filterSNPs(SNPset = d15.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d15.2.filt <- runQTLseqAnalysis(SNPset = d15.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d15.2.filt <- runGprimeAnalysis(d15.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d15.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d15.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day5.50nM-control.rep2.csv",row.names=FALSE)


d15.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0225",lowBulk = "AB_BC_0229")
d15.3.filt <-filterSNPs(SNPset = d15.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d15.3.filt <- runQTLseqAnalysis(SNPset = d15.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d15.3.filt <- runGprimeAnalysis(d15.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d15.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d15.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day5.100nM-control.rep1.csv",row.names=FALSE)
d15.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0226",lowBulk = "AB_BC_0230")
d15.4.filt <-filterSNPs(SNPset = d15.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d15.4.filt <- runQTLseqAnalysis(SNPset = d15.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d15.4.filt <- runGprimeAnalysis(d15.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d15.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d15.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day5.100nM-control.rep2.csv",row.names=FALSE)

d15.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0225",lowBulk = "AB_BC_0227")
d15.5.filt <-filterSNPs(SNPset = d15.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d15.5.filt <- runQTLseqAnalysis(SNPset = d15.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d15.5.filt <- runGprimeAnalysis(d15.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d15.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d15.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day5.250nM-control.rep1.csv",row.names=FALSE)
d15.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0226",lowBulk = "AB_BC_0228")
d15.6.filt <-filterSNPs(SNPset = d15.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d15.6.filt <- runQTLseqAnalysis(SNPset = d15.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d15.6.filt <- runGprimeAnalysis(d15.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d15.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d15.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day5.250nM-control.rep2.csv",row.names=FALSE)


################################################################################################################
################################################################################################################

# M26 96h day10 - control
d16.1 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0529",lowBulk = "AB_BC_0535")
d16.1.filt <-filterSNPs(SNPset = d16.1,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d16.1.filt <- runQTLseqAnalysis(SNPset = d16.1.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d16.1.filt <- runGprimeAnalysis(d16.1.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d16.1.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d16.1.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day10.50nM-control.rep1.csv",row.names=FALSE)
d16.2 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0530",lowBulk = "AB_BC_0536")
d16.2.filt <-filterSNPs(SNPset = d16.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d16.2.filt <- runQTLseqAnalysis(SNPset = d16.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d16.2.filt <- runGprimeAnalysis(d16.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d16.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d16.2.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day10.50nM-control.rep2.csv",row.names=FALSE)


d16.3 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0529",lowBulk = "AB_BC_0533")
d16.3.filt <-filterSNPs(SNPset = d16.3,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d16.3.filt <- runQTLseqAnalysis(SNPset = d16.3.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d16.3.filt <- runGprimeAnalysis(d16.3.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d16.3.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d16.3.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day10.100nM-control.rep1.csv",row.names=FALSE)
d16.4 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0530",lowBulk = "AB_BC_0534")
d16.4.filt <-filterSNPs(SNPset = d16.4,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d16.4.filt <- runQTLseqAnalysis(SNPset = d16.4.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d16.4.filt <- runGprimeAnalysis(d16.4.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d16.4.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d16.4.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day10.100nM-control.rep2.csv",row.names=FALSE)

d16.5 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0529",lowBulk = "AB_BC_0531")
d16.5.filt <-filterSNPs(SNPset = d16.5,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d16.5.filt <- runQTLseqAnalysis(SNPset = d16.5.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d16.5.filt <- runGprimeAnalysis(d16.5.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d16.5.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d16.5.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day10.250nM-control.rep1.csv",row.names=FALSE)
d16.6 <-importFromGATK(file = "BSA5.NF54-NHP4026.SNP.calledINparents.coreGenome2.table",highBulk = "AB_BC_0530",lowBulk = "AB_BC_0532")
d16.6.filt <-filterSNPs(SNPset = d16.6,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d16.6.filt <- runQTLseqAnalysis(SNPset = d16.6.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))
d16.6.filt <- runGprimeAnalysis(d16.6.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = d16.6.filt, var = "Gprime", plotThreshold = FALSE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()+ylim(0,100)
write.csv(d16.6.filt[,c(1:4,1441:1465)], file = "QTLs.CQ/M26.96h.day10.250nM-control.rep2.csv",row.names=FALSE)




