##This file contains all the scripts used for RNA expression, copy number and clinical data analysis 

###
ucecRNA = read.table("./gdac.broadinstitute.org_UCEC.mRNAseq_Preprocess.Level_3.2016012800.0.0/UCEC.uncv2.mRNAseq_RSEM_normalized_log2.txt",sep = "\t",header = T, row.names = 1)
dim(ucecRNA) #18536   581
ovRNA = read.table("./OV/gdac.broadinstitute.org_OV.mRNAseq_Preprocess.Level_3.2016012800.0.0/OV.uncv2.mRNAseq_RSEM_normalized_log2.txt", sep = "\t", row.names = 1, header = T)
dim(ovRNA)
range(ovRNA,na.rm = T)
commonGenes = intersect(rownames(ucecRNA),rownames(ovRNA))
commonGenes = commonGenes[!grepl(commonGenes,pattern = "\\?")]
ucecRNA = ucecRNA[commonGenes, ]
ovRNA = ovRNA[commonGenes,]
range(ucecRNA,na.rm = T)
range(range(ucecRNA,na.rm = T))
table(complete.cases(ucecRNA))
table(complete.cases(ovRNA))
library(DMwR)
saveRDS(ucecRNA, file = "ucecRNA.rds")
saveRDS(ovRNA, file = "ovRNA.rds")
# ucecRNA = knnImputation(ucecRNA)
# ovRNA = knnImputation(ovRNA)
ucecRNA = readRDS("./afterKNNucecRNA.rds")
ovRNA = readRDS("./afterKNNovRNA.rds")
# table(substr(colnames(ucecRNA),13,16))
# table(substr(colnames(ovRNA),13,16))
colnames(ovRNA) = gsub(colnames(ovRNA), pattern = "\\.", replacement = "-")

ovRNA = ovRNA[,grepl(colnames(ovRNA),pattern = "01$")]
colnames(ovRNA) = gsub(colnames(ovRNA), pattern = "-01", replacement = "")

ucecHis = setNames(suhasUCEC$histological_type,rownames(suhasUCEC))
table(ucecHis) #now from suhas

table(substr(colnames(ucecRNA), 13,15))
table(substr(colnames(ovRNA), 13,15))

ucecRNA = ucecRNA[,grepl(colnames(ucecRNA), pattern = "01$")]
colnames(ucecRNA) = substr(colnames(ucecRNA), 1, 12)
colnames(ucecRNA) = gsub(colnames(ucecRNA), pattern = "\\.", replacement = "-")
table(names(ucecHis) %in% colnames(ucecRNA)) #TRUE 523
ucecRNA = ucecRNA[,intersect(names(ucecHis), colnames(ucecRNA))]
ucecHis = ucecHis[intersect(names(ucecHis), colnames(ucecRNA))]

totalRNA = cbind(ucecRNA,ovRNA)
totalRNA = as.data.frame(t(totalRNA))

totalRNA$anno = c(ucecHis, rep("OV", ncol(ovRNA)))
library(ggfortify)
autoplot(prcomp(totalRNA[,1:(ncol(totalRNA)-1)]), data = totalRNA, colour = "anno", frame = TRUE, frame.type = 'norm')
autoplot(prcomp(totalRNA[,1:(ncol(totalRNA)-1)]), data = totalRNA, colour = "anno")

####performing sample-wise spearman correlation for all possible sample pairs
selectGenes = apply(totalRNA[,-ncol(totalRNA)], 2, var)
selectGenes = sort(selectGenes, decreasing = TRUE)
selectGenes = names(selectGenes)[1:1000]

endoSerousCor = cor(t(totalRNA[totalRNA$anno=="Endometrioid",selectGenes]), t(totalRNA[totalRNA$anno=="Serous",selectGenes]), method = "spearman")
mean(endoSerousCor)
sd(endoSerousCor)

endoOVCor = cor(t(totalRNA[totalRNA$anno=="Endometrioid",selectGenes]), t(totalRNA[totalRNA$anno=="OV",selectGenes]), method = "spearman")
mean(endoOVCor)
sd(endoOVCor)

SerousOVCor = cor(t(totalRNA[totalRNA$anno=="Serous",selectGenes]), t(totalRNA[totalRNA$anno=="OV",selectGenes]), method = "spearman")
mean(SerousOVCor)
sd(SerousOVCor)

t.test(endoSerousCor,SerousOVCor)
t.test(endoSerousCor,endoOVCor)
t.test(SerousOVCor,endoOVCor)

endoCor = cor(t(totalRNA[totalRNA$anno=="Endometrioid",selectGenes]), t(totalRNA[totalRNA$anno=="Endometrioid",selectGenes]), method = "spearman")
mean(endoCor)
sd(endoCor)

OVCor = cor(t(totalRNA[totalRNA$anno=="OV",selectGenes]), t(totalRNA[totalRNA$anno=="OV",selectGenes]), method = "spearman")
mean(OVCor)
sd(OVCor)

SerousCor = cor(t(totalRNA[totalRNA$anno=="Serous",selectGenes]), t(totalRNA[totalRNA$anno=="Serous",selectGenes]), method = "spearman")
mean(SerousCor)
sd(SerousCor)

sampleWiseCorDf = data.frame(cor = c(endoCor,SerousCor,OVCor,endoSerousCor,SerousOVCor,endoOVCor), 
           type = rep(c("Endometroid EC","Serous-like EC", "Serous OvCa","Endometroid EC vs. Serous-like EC","Serous-like EC vs. Serous OvCa","Endometroid EC vs. Serous OvCa"),
               c(length(endoCor),length(SerousCor), length(OVCor), length(endoSerousCor), length(SerousOVCor), length(endoOVCor))))
sampleWiseCorDf$type = factor(sampleWiseCorDf$type, levels = unique(sampleWiseCorDf$type))
ggplot(sampleWiseCorDf, aes(type, cor))+geom_boxplot()
####check the expression of tissue-specific gene
endoSpGene = read.table("./uterus_genes.txt", sep = "\t", header = T, stringsAsFactors = F, skip = 1)
ovSpGene = read.table("ovary_genes.txt", sep = "\t", header = T, stringsAsFactors = F,skip = 1)
tissueGenes = list(uterus = endoSpGene$Gene_Symbol, ovary = ovSpGene$Gene_Symbol)
totalRNA = t(totalRNA)
rownames(totalRNA) = gsub(rownames(totalRNA), pattern = "\\|.*", replacement = "")
library(GSVA)
tissueGeneScore = gsva(totalRNA, tissueGenes, "ssgsea")
tissueGeneScore = as.data.frame(t(tissueGeneScore))
tissueGeneScore = tissueGeneScore[c(endosurvsamples,seroussurvsamples,ovsurvsamples),]
tissueGeneScore$group = rep(c("Endometroid EC","Serous-like EC","Serous OV"),
                            c(length(endosurvsamples), length(seroussurvsamples), length(ovsurvsamples)))
tissueGeneScore$group = factor(tissueGeneScore$group, levels = c("Endometroid EC","Serous-like EC","Serous OV"))
ggplot(tissueGeneScore)+geom_boxplot(aes(x=group,y=uterus))+ylab("Uterus-specific gene expression")
ggplot(tissueGeneScore)+geom_boxplot(aes(x=group,y=ovary))+ylab("Ovary-specific gene expression")
t.test(tissueGeneScore[tissueGeneScore$group == "Endometroid EC","uterus"],
       tissueGeneScore[tissueGeneScore$group == "Serous-like EC","uterus"])
t.test(tissueGeneScore[tissueGeneScore$group == "Endometroid EC","uterus"],
       tissueGeneScore[tissueGeneScore$group == "Serous OV","uterus"])

t.test(tissueGeneScore[tissueGeneScore$group == "Endometroid EC","ovary"],
       tissueGeneScore[tissueGeneScore$group == "Serous-like EC","ovary"])
t.test(tissueGeneScore[tissueGeneScore$group == "Endometroid EC","ovary"],
       tissueGeneScore[tissueGeneScore$group == "Serous OV","ovary"])
###

ovSCNA$sample = substr(ovSCNA$sample, 1, 12)
ucecSCNA$sample = substr(ucecSCNA$sample,1,12)

ovSCNA = ovSCNA[ovSCNA$sample %in% colnames(ovRNA),]
length(intersect(ovSCNA$sample,colnames(ovRNA))) #289
length(intersect(ucecSCNA$sample,colnames(ucecRNA))) #345

ucecSCNA = ucecSCNA[ucecSCNA$sample %in% colnames(ucecRNA),]
source("~/OneDrive/deconv/scrips/all_in_one_script.r")

plotCNV(ucecSCNA[ucecSCNA$sample %in% names(ucecHis)[ucecHis == "Endometrioid"],])
plotCNV(ucecSCNA[ucecSCNA$sample %in% names(ucecHis)[ucecHis == "Serous"],])
plotCNV(ovSCNA)

length(unique(ucecSCNA$sample[ucecSCNA$sample %in% names(ucecHis)[ucecHis == "Endometrioid"]]))
length(unique(ucecSCNA$sample[ucecSCNA$sample %in% names(ucecHis)[ucecHis == "Serous"]]))
length(unique(ovSCNA$sample))

delDiff = list()
ampDiff = list()
colnames(ucecCNVDel) = substr(colnames(ucecCNVDel), 1, 12)
colnames(ovCNVDel) = substr(colnames(ovCNVDel),1,12)
colnames(ucecCNVAmp) = substr(colnames(ucecCNVAmp), 1, 12)
colnames(ovCNVDel) = substr(colnames(ovCNVDel), 1, 12)

ucecEnSamples = unique(ucecSCNA$sample[ucecSCNA$sample %in% names(ucecHis)[ucecHis == "Endometrioid"]])
ucecSerousSamples = unique(ucecSCNA$sample[ucecSCNA$sample %in% names(ucecHis)[ucecHis == "Serous"]])
ovSamples = unique(ovSCNA$sample)
for(arm in rownames(ucecCNVDel)){
  delDiff[[arm]] = t.test(as.numeric(ucecCNVDel[arm, ucecSerousSamples]),as.numeric(ovCNVDel[arm, ]))
}
delDiff = sapply(delDiff, function(x) x$p.value) 
which(delDiff < 0.01)
delDiffArm = c("chr6q_del","chr13q_del","chr7p_del","chr18q_del")
plotDel = cbind(ucecCNVDel[delDiffArm,ucecSerousSamples],ovCNVDel[delDiffArm,])
plotDel = as.data.frame(t(plotDel))
plotDel$type = c(rep("UCEC",length(ucecSerousSamples)),rep("OV",ncol(ovCNVDel)))
plotDel$type = factor(plotDel$type, levels = c("UCEC","OV"))
plotDel = plotDel %>% gather(key = "Arm", value = "log2Ratio", 1:4)
ggplot(plotDel,aes(Arm, log2Ratio, fill = type))+geom_boxplot()

for(arm in rownames(ucecCNVAmp)){
  ampDiff[[arm]] = t.test(as.numeric(ucecCNVAmp[arm, ucecSerousSamples]),as.numeric(ovCNVAmp[arm, ]))
}
ampDiff = sapply(ampDiff, function(x) x$p.value) 
which(ampDiff < 0.01)
AmpDiffArm = c("chr1q_amp","chr11q_amp","chr12p_amp")
plotAmp = cbind(ucecCNVAmp[AmpDiffArm,ucecSerousSamples],ovCNVAmp[AmpDiffArm,])
plotAmp = as.data.frame(t(plotAmp))
plotAmp$type = c(rep("UCEC",length(ucecSerousSamples)),rep("OV",ncol(ovCNVDel)))
plotAmp$type = factor(plotAmp$type, levels = c("UCEC","OV"))
plotAmp = plotAmp %>% gather(key = "Arm", value = "log2Ratio", 1:3)
ggplot(plotAmp,aes(Arm, log2Ratio, fill = type))+geom_boxplot()

ucecGeneLevel = fread("./gdac.broadinstitute.org_UCEC-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_data_by_genes.txt",
                           sep = "\t")
ovGeneLevel = fread("./gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_data_by_genes.txt",
                      sep = "\t")
ucecGeneLevelB = fread("./gdac.broadinstitute.org_UCEC-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt",
                       sep = "\t")
ovGeneLevelB = fread("./gdac.broadinstitute.org_OV-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0/all_thresholded.by_genes.txt",
                       sep = "\t")
#only keep the chromosome 4
ucecGeneLevel = as.data.frame(ucecGeneLevel)
ovGeneLevel = as.data.frame(ovGeneLevel)
ucecGeneLevelB = as.data.frame(ucecGeneLevelB)
ovGeneLevelB = as.data.frame(ovGeneLevelB)
#ovGeneLevelB = as.data.frame(ovGeneLevelB[grepl(Cytoband, pattern = "^4")])
rownames(ucecGeneLevel) = ucecGeneLevel$`Gene Symbol`
rownames(ucecGeneLevelB) = ucecGeneLevelB$`Gene Symbol`
rownames(ovGeneLevel) = ovGeneLevel$`Gene Symbol`
rownames(ovGeneLevelB) = ovGeneLevelB$`Gene Symbol`

ucecGeneLevel = ucecGeneLevel[grepl(colnames(ucecGeneLevel), pattern = "-01A")]
ucecGeneLevelB = ucecGeneLevelB[grepl(colnames(ucecGeneLevelB), pattern = "-01A")]
ovGeneLevel = ovGeneLevel[grepl(colnames(ovGeneLevel), pattern = "-01A")]
ovGeneLevelB = ovGeneLevelB[grepl(colnames(ovGeneLevelB), pattern = "-01A")]
ovGeneLevelB = ovGeneLevelB[grepl(colnames(ovGeneLevelB), pattern = "-01A")]

colnames(ucecGeneLevel) = substr(colnames(ucecGeneLevel), 1, 12)
colnames(ucecGeneLevelB) = substr(colnames(ucecGeneLevelB), 1, 12)
colnames(ovGeneLevel) = substr(colnames(ovGeneLevel), 1, 12)
colnames(ovGeneLevelB) = substr(colnames(ovGeneLevelB), 1, 12)

endoGeneLevel= ucecGeneLevel[, rownames(totalCNV)[totalCNV$anno == "Endometrioid"]]
serousGeneLevel= ucecGeneLevel[, rownames(totalCNV)[totalCNV$anno == "Serous"]]
ovGeneLevel= ovGeneLevel[, rownames(totalCNV)[totalCNV$anno == "OV"]]
endoGeneLevelB= ucecGeneLevelB[, rownames(totalCNV)[totalCNV$anno == "Endometrioid"]]
serousGeneLevelB= ucecGeneLevelB[, rownames(totalCNV)[totalCNV$anno == "Serous"]]
ovGeneLevelB= ovGeneLevelB[, rownames(totalCNV)[totalCNV$anno == "OV"]]

all(rownames(endoGeneLevel) == rownames(serousGeneLevel))
all(rownames(endoGeneLevel) == rownames(ovGeneLevelB))
chr4CN_diff = lapply(1:nrow(endoGeneLevel) , function(x)t.test(as.numeric(serousGeneLevel[x,]),as.numeric(ovGeneLevelB[x,])))
names(chr4CN_diff) = rownames(endoGeneLevel)
chr4CN_diff_p = sapply(chr4CN_diff, function(x)x$p.value)
chr4CN_diff_p = sort(chr4CN_diff_p)
chr4CN_diff_p[1:50]
cosmic = fread("~/OneDrive/Census_allWed Sep 18 04_11_25 2019.tsv")
cosmicSuppressor = cosmic[grepl(`Role in Cancer`, pattern = "TSG"),`Gene Symbol`]
chr4CN_diff_p_sig=names(chr4CN_diff_p[intersect(cosmicSuppressor,names(chr4CN_diff))])
chr4qPlot = cbind(endoGeneLevel[chr4CN_diff_p_sig,],
                  serousGeneLevel[chr4CN_diff_p_sig,],
                  ovGeneLevel[chr4CN_diff_p_sig,])
chr4qPlot = t(chr4qPlot)
chr4qPlot = as.data.frame(chr4qPlot)
chr4qPlot$anno = rep(c("Endo","Serous","OV"), c(ncol(endoGeneLevel), ncol(serousGeneLevel), ncol(ovGeneLevel)))
chr4qPlot$anno = factor(chr4qPlot$anno, levels = c("Endo","Serous","OV"))
chr4qPlot = chr4qPlot %>% gather(key = "gene", value = "cn", 1:(ncol(chr4qPlot)-1))
ggplot(chr4qPlot, aes(gene, cn, col = anno))+
  geom_boxplot(outlier.size = 0.5)+
  theme(axis.text.x = element_text(angle = -90))

endoChr4Del = apply(endoGeneLevelB[chr4CN_diff_p_sig,], 1, function(x)sum(x< 0)/length(x))
SerousChr4Del = apply(serousGeneLevelB[chr4CN_diff_p_sig,], 1, function(x)sum(x< 0)/length(x))
OvChr4Del = apply(ovGeneLevelB[chr4CN_diff_p_sig,], 1, function(x)sum(x < 0)/length(x))
                  
chr4qPlot2 = data.frame(gene = rep(chr4CN_diff_p_sig,3), 
                        proportion = c(endoChr4Del,SerousChr4Del,OvChr4Del),
                        anno = rep(c("Endo","Serous","OV"), each=length(chr4CN_diff_p_sig)))
chr4qPlot2$anno = factor(chr4qPlot2$anno, levels = c("Endo","Serous","OV"))

ggplot(chr4qPlot2, aes(gene, proportion, fill = anno))+
  geom_bar(stat = "identity", position=position_dodge())+
  theme(axis.text.x = element_text(angle = -90))


#####
library(TCGAbiolinks)
dim(ucecHi) 
colData(ucecHi)$subtype_histology 
dim(ucecGA)

#the ucec clinic info from linkedOmics
suhasUCEC = read.table("./Human__TCGA_UCEC__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi",
                       sep = "\t", row.names = 1, header = T, stringsAsFactors = F)
suhasUCEC = as.data.frame(t(suhasUCEC))
table(suhasUCEC$histological_type)
dim(suhasUCEC)
suhasUCEC = suhasUCEC[suhasUCEC$histological_type!= "mixedserousandendometrioid",]
suhasUCEC$histological_type = gsub(suhasUCEC$histological_type, pattern = "endometrioidendometrialadenocarcinoma", replacement = "Endometrioid")
suhasUCEC$histological_type = gsub(suhasUCEC$histological_type, pattern = "serousendometrialadenocarcinoma", replacement = "Serous")
rownames(suhasUCEC) = gsub(rownames(suhasUCEC), pattern = "\\.", replacement = "-")

####
ucecLinc = read.table("~/ucecTCGA/ucec_clinic.txt", sep = "\t", header = TRUE, stringsAsFactors = F)
ucecLinc = ucecLinc[!duplicated(ucecLinc$bcr_patient_barcode),]
dim(ucecLinc)
ucecLinc$stage_event_clinical_stage
table(ucecLinc$histological_type)
SerousClinc = ucecLinc[ucecLinc$histological_type == "Serous endometrial adenocarcinoma",]
endoClinc = ucecLinc[ucecLinc$histological_type == "Endometrioid endometrial adenocarcinoma",]
table(SerousClinc$stage_event_clinical_stage)
SerousClinc$stage = ifelse(grepl(SerousClinc$stage_event_clinical_stage,pattern = "(III)|(IV)"), "late","early")
endoClinc$stage = ifelse(grepl(endoClinc$stage_event_clinical_stage,pattern = "(III)|(IV)"), "late","early")
table(endoClinc$stage) #326    85 
table(SerousClinc$stage) #54    61 

query <- GDCquery(project = "TCGA-OV",
                  data.category = "Clinical",
                  file.type = "xml")
GDCdownload(query)
ovClinic <- GDCprepare_clinic(query, clinical.info = "patient")
dim(ovClinic) #587  59
ovClinic = ovClinic[!duplicated(ovClinic$bcr_patient_barcode),]
table(ovClinic$stage_event_clinical_stage)
ovClinic$stage = ifelse(grepl(ovClinic$stage_event_clinical_stage, pattern = "(III)|(IV)"),"late","early")
table(ovClinic$stage)
ovSampleSub = intersect(colnames(ovRNA), ovClinic$bcr_patient_barcode)
length(ovSampleSub) #303
ovClinicSub = ovClinic[ovClinic$bcr_patient_barcode %in% ovSampleSub,]
table(ovClinicSub$stage) #24 279


