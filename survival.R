#this file includes the scripts for survival mutation analysis

table(endoClinc$bcr_patient_barcode %in% rownames(suhasUCEC)) #all TRUE
#use suhasUCEC for survival
table(suhasUCEC$histological_type)
suhasUCEC$overall_survival = as.numeric(suhasUCEC$overall_survival)
suhasUCEC$status = as.numeric(suhasUCEC$status)
suhasUCEC = suhasUCEC[!is.na(suhasUCEC$overall_survival),]
suhasEndo = suhasUCEC[suhasUCEC$histological_type == "Endometrioid",]
suhasSerous = suhasUCEC[suhasUCEC$histological_type == "Serous",]

suhasOV = read.table("./Human__TCGA_OV__MS__Clinical__Clinical__01_28_2016__BI__Clinical__Firehose.tsi.txt",
                     sep = "\t", header = TRUE, stringsAsFactors = F, row.names = 1 )
suhasOV = t(suhasOV)
suhasOV = as.data.frame(suhasOV)
rownames(suhasOV) = gsub(rownames(suhasOV), pattern = "\\.", replacement = "-")
suhasOV$overall_survival = as.numeric(suhasOV$overall_survival)
suhasOV$status = as.numeric(suhasOV$status)
dim(suhasOV)
suhasOV = suhasOV[!is.na(suhasOV$overall_survival),]

library(survcomp)
ovRNA = ovRNA[, grepl(colnames(ovRNA), pattern = "01$")]
colnames(ovRNA) = substr(colnames(ovRNA), 1, 12)
colnames(ovRNA) = gsub(colnames(ovRNA), pattern = "\\.", replacement = "-")
ovsurvsamples = intersect(rownames(suhasOV), colnames(ovRNA))
length(ovsurvsamples) #294

endosurvsamples = intersect(rownames(suhasEndo), colnames(ucecRNA))
length(endosurvsamples) #395
seroussurvsamples = intersect(rownames(suhasSerous), colnames(ucecRNA))
length(seroussurvsamples) #110

ovSurv = rep(NA, nrow(ovRNA))
names(ovSurv) = rownames(ovRNA)
endoSurv = rep(NA, nrow(ucecRNA))
names(endoSurv) = rownames(ucecRNA)
serousSurv = endoSurv


for(gene in rownames(ovRNA))
  ovSurv[gene] = concordance.index(x = as.numeric(ovRNA[gene,ovsurvsamples]), 
                                   surv.time = suhasOV[ovsurvsamples,"overall_survival"],
                                   surv.event = suhasOV[ovsurvsamples,"status"])
ovSurv = unlist(ovSurv)
ovCox = list()

for(gene in rownames(ucecRNA)){
  endoSurv[gene] = concordance.index(x = as.numeric(ucecRNA[gene,endosurvsamples]), 
                                   surv.time = suhasEndo[endosurvsamples,"overall_survival"],
                                   surv.event = suhasEndo[endosurvsamples,"status"])
  serousSurv[gene] = concordance.index(x = as.numeric(ucecRNA[gene,seroussurvsamples]), 
                                     surv.time = suhasSerous[seroussurvsamples,"overall_survival"],
                                     surv.event = suhasSerous[seroussurvsamples,"status"])
}
serousSurv = unlist(serousSurv)
endoSurv = unlist(endoSurv)

ovSurv = ovSurv[order(abs(ovSurv-0.5), decreasing = TRUE)]
serousSurv = serousSurv[order(abs(serousSurv-0.5), decreasing = TRUE)]
endoSurv = endoSurv[order(abs(endoSurv-0.5), decreasing = TRUE)]
intersect(names(ovSurv)[ovSurv<0.4],names(serousSurv)[serousSurv<0.4])
intersect(names(endoSurv)[endoSurv>0.6],names(serousSurv)[serousSurv>0.6])
intersect(names(endoSurv)[1:100],names(serousSurv)[1:100])
intersect(names(endoSurv)[1:100],names(ovSurv)[1:100])

source('~/OneDrive/deconv/scrips/survivalByMedian.r')
survCombined = Surv(time = c(suhasEndo$overall_survival, suhasSerous$overall_survival, suhasOV$overall_survival),
                      event = c(suhasEndo$status, suhasSerous$status, suhasOV$status))
names(survCombined) = c(rownames(suhasEndo), rownames(suhasSerous), rownames(suhasOV))
survCombinedDf = data.frame(surv=survCombined, group = rep(c("Endo","Serous","OV"),c(nrow(suhasEndo), nrow(suhasSerous), nrow(suhasOV))))
survCombinedDf$caseID = c(rownames(suhasEndo), rownames(suhasSerous), rownames(suhasOV))
survCombinedDf$group = factor(survCombinedDf$group, levels = c("Endo","Serous","OV"))
subtypeFit <- survfit(survCombined ~ group, data = survCombinedDf)
library(survminer)
ggsurvplot(subtypeFit, pval = TRUE)

survCombinedList = split(survCombinedDf,f = survCombinedDf$group)
ovSurvCox = list()
for(gene in rownames(ovRNA)){
  ovSurvCox[[gene]] = coxph(survCombinedList$OV[ovsurvsamples,"surv"] ~ as.numeric(ovRNA[gene,ovsurvsamples]))
}
ovSurvCoxP = sapply(ovSurvCox, function(x)summary(x)$logtest[3])
ovSurvCoxP = sort(ovSurvCoxP)

endoSurvCox = list()
serousSurCox = list()
for(gene in rownames(ucecRNA)){
  serousSurCox[[gene]] = coxph(survCombinedList$Serous[seroussurvsamples,"surv"] ~ as.numeric(ucecRNA[gene,seroussurvsamples]))
  endoSurvCox[[gene]] = coxph(survCombinedList$Endo[endosurvsamples,"surv"] ~ as.numeric(ucecRNA[gene,endosurvsamples]))
}

serousSurCoxP = sapply(serousSurCox, function(x)summary(x)$logtest[3])
serousSurCoxP = sort(serousSurCoxP)
endoSurvCoxP = sapply(endoSurvCox, function(x)summary(x)$logtest[3])
endoSurvCoxP = sort(endoSurvCoxP)
twoEC_commonSurvGenes = intersect(names(serousSurCoxP)[serousSurCoxP<0.05], names(endoSurvCoxP)[endoSurvCoxP<0.05])
intersect(names(serousSurCoxP)[serousSurCoxP<0.05], names(ovSurvCoxP)[ovSurv<0.05])
twoEC_commonSurvGenes_logRank_endo = sapply(twoEC_commonSurvGenes, function(x)survByMedian(scoreVec = setNames(as.numeric(ucecRNA[x,endosurvsamples]), endosurvsamples),
                                                                               surv = data.frame(time = suhasEndo[endosurvsamples,"overall_survival"], 
                                                                                                 event = suhasEndo[endosurvsamples,"status"], 
                                                                                                 row.names = endosurvsamples), 
                                                                               plot = F))
twoEC_commonSurvGenes_logRank_serous = sapply(twoEC_commonSurvGenes, function(x)survByMedian(scoreVec = setNames(as.numeric(ucecRNA[x,seroussurvsamples]), seroussurvsamples),
                                                                                           surv = data.frame(time = suhasSerous[seroussurvsamples,"overall_survival"], 
                                                                                                             event = suhasSerous[seroussurvsamples,"status"], 
                                                                                                             row.names = seroussurvsamples), 
                                                                                           plot = F))
twoEC_commonSurvGene_cIdx = sapply(twoEC_commonSurvGenes, function(x)survByMedian(scoreVec = setNames(as.numeric(ucecRNA[x,seroussurvsamples]), seroussurvsamples),
                                                                                  surv = data.frame(time = suhasSerous[seroussurvsamples,"overall_survival"], 
                                                                                                    event = suhasSerous[seroussurvsamples,"status"], 
                                                                                                    row.names = seroussurvsamples), 
                                                                                  plot = F)$himsc.cIndex)
twoEC_commonSurvGenes_df = data.frame(symbol = gsub(twoEC_commonSurvGenes, pattern = "\\|.*", replacement = ""),
                                       entrezID = gsub(twoEC_commonSurvGenes, pattern = ".*\\|", replacement = ""),
                                        CoxPH_p_value_endo = endoSurvCoxP[twoEC_commonSurvGenes],
                                       CoxPH_p_value_serous = serousSurCoxP[twoEC_commonSurvGenes],
                                       log_rank_pvalue_endo = twoEC_commonSurvGenes_logRank_endo,
                                      log_rank_pvalue_serous = twoEC_commonSurvGenes_logRank_serous,
                                      prog_trend = ifelse(twoEC_commonSurvGene_cIdx<0.5, "good","bad"))
data.table::fwrite(twoEC_commonSurvGenes_df, file = "commonSurvival_genes_table_2.txt", sep = "\t")

survByMedian(scoreVec = setNames(as.numeric(ucecRNA["PHKA1",endosurvsamples]), endosurvsamples),
             surv = data.frame(time = suhasEndo[endosurvsamples,"overall_survival"], 
                               event = suhasEndo[endosurvsamples,"status"], row.names = endosurvsamples),
             plot = T)
survByMedian(scoreVec = setNames(as.numeric(ucecRNA["PHKA1", seroussurvsamples]), seroussurvsamples),
             surv = data.frame(time = suhasSerous[seroussurvsamples,"overall_survival"], 
                               event = suhasSerous[seroussurvsamples,"status"], row.names = seroussurvsamples),
             plot = T)
survByMedian(scoreVec = setNames(as.numeric(ovRNA["PHKA1", ovsurvsamples]), ovsurvsamples),
             surv = data.frame(time = suhasOV[ovsurvsamples,"overall_survival"], 
                               event = suhasOV[ovsurvsamples,"status"], row.names = ovsurvsamples),
             plot = T)

survByMedian(scoreVec = setNames(as.numeric(ucecRNA["CXCR5",endosurvsamples]), endosurvsamples),
             surv = data.frame(time = suhasEndo[endosurvsamples,"overall_survival"], 
                               event = suhasEndo[endosurvsamples,"status"], row.names = endosurvsamples),
             plot = T)
survByMedian(scoreVec = setNames(as.numeric(ucecRNA["CXCR5", seroussurvsamples]), seroussurvsamples),
             surv = data.frame(time = suhasSerous[seroussurvsamples,"overall_survival"], 
                               event = suhasSerous[seroussurvsamples,"status"], row.names = seroussurvsamples),
             plot = T)
survByMedian(scoreVec = setNames(as.numeric(ovRNA["CXCR5", ovsurvsamples]), ovsurvsamples),
             surv = data.frame(time = suhasOV[ovsurvsamples,"overall_survival"], 
                               event = suhasOV[ovsurvsamples,"status"], row.names = ovsurvsamples),
             plot = T)




ovCoxCoef = sapply(ovSurvCox, function(x)summary(x)$coef[2])
ovCoxCoef = sort(ovCoxCoef)
ovCoxCoefDf = data.frame(gene = gsub(names(ovCoxCoef), pattern = "\\|.*", replacement = ""),
                         coef = ovCoxCoef)
write.table(ovCoxCoefDf, file = "ovCoxRank.rnk", sep = "\t", row.names = F, col.names = F, quote = F)

endoCoxCoef = sapply(endoSurvCox, function(x)summary(x)$coef[2])
endoCoxCoef = sort(endoCoxCoef)
endoCoxCoefDf = data.frame(gene = gsub(names(endoCoxCoef), pattern = "\\|.*", replacement = ""),
                         coef = endoCoxCoef)
write.table(endoCoxCoefDf, file = "endoCoxRank.rnk", sep = "\t", row.names = F, col.names = F, quote = F)

serousCoxCoef = sapply(serousSurCox, function(x)summary(x)$coef[2])
serousCoxCoef = sort(serousCoxCoef)
serousCoxCoefDf = data.frame(gene = gsub(names(serousCoxCoef), pattern = "\\|.*", replacement = ""),
                         coef = serousCoxCoef)
write.table(serousCoxCoefDf, file = "serousCoxRank.rnk", sep = "\t", row.names = F, col.names = F, quote = F)

ovRNA = ovRNA[!duplicated(gsub(pattern = "\\|.*", replacement = "", rownames(ovRNA))),]
rownames(ovRNA) = gsub(pattern = "\\|.*", replacement = "", rownames(ovRNA))
ucecRNA = ucecRNA[!duplicated(gsub(pattern = "\\|.*", replacement = "", rownames(ucecRNA))),]
rownames(ucecRNA) = gsub(pattern = "\\|.*", replacement = "", rownames(ucecRNA))

###output survival data
table(endoSurvCoxP<0.01) #623
table(serousSurCoxP<0.01) #101
table(ovSurvCoxP<0.01) #255
names(endoSurvCoxP) = gsub(names(endoSqurvCoxP), pattern = "\\.pvalue", replacement = "")
endoCoxCoefDf$p.value = endoSurvCoxP[rownames(endoCoxCoefDf)]

names(serousSurCoxP) = gsub(names(serousSurCoxP), pattern = "\\.pvalue", replacement = "")
serousCoxCoefDf$p.value = serousSurCoxP[rownames(serousCoxCoefDf)]

names(ovSurvCoxP) = gsub(names(ovSurvCoxP), pattern = "\\.pvalue", replacement = "")
ovCoxCoefDf$p.value = ovSurvCoxP[rownames(ovCoxCoefDf)]

endoCoxCoefDf = endoCoxCoefDf[order(endoCoxCoefDf$p.value),]
serousCoxCoefDf = serousCoxCoefDf[order(serousCoxCoefDf$p.value),]
ovCoxCoefDf = ovCoxCoefDf[order(ovCoxCoefDf$p.value),]
combinedCoxCoefDf = rbind(endoCoxCoefDf[1:623,],serousCoxCoefDf[1:101,], ovCoxCoefDf[1:255,])
combinedCoxCoefDf$cancer = rep(c("Endometroid EC","Serous-like EC","serous OvCa"),c(623, 101, 255))
write.table(combinedCoxCoefDf, sep = "\t", row.names = F, col.names = T, quote = F, file = "Cox_Survival_gene.txt")

####
write.table(ovRNA, file = "ovRNA4ciber.txt", sep = "\t", col.names = T, row.names = T, quote = F)
write.table(ucecRNA, file = "ucecRNA4ciber.txt", sep = "\t", col.names = T, row.names = T, quote = F)

