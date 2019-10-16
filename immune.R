#This file contains all scripts for immune-related analysis
library(estimate)
ovEst = estimateScore(input.ds = "./ovRNA4ciber.gct.txt", output.ds = "./ovEst.gct")
ovEst = CePa::read.gct("./ovEst.gct")
ovEst = t(ovEst)
ovEst = as.data.frame(ovEst)
ucecEst = estimateScore(input.ds = "./ucecRNA4ciber.gct.txt", output.ds = "./ucecEst.gct")
ucecEst = CePa::read.gct("./ucecEst.gct")
ucecEst = t(ucecEst)
ucecEst = as.data.frame(ucecEst)
rownames(ucecEst) = gsub(rownames(ucecEst), pattern = "\\.", replacement = "-")
rownames(ovEst) = gsub(rownames(ovEst), pattern = "\\.", replacement = "-")
#table(rownames(ucecEst) %in% colnames(ucecRNA))
boxplot(ucecEst$StromalScore[rownames(ucecEst) %in% rownames(suhasEndo)],
        ucecEst$StromalScore[rownames(ucecEst) %in% rownames(suhasSerous)],
        ovEst$StromalScore, names = c("Endometroid EC", "Serous-like EC", "Serous OvCa"), 
        ylab = "Stromal score", las=2)

boxplot(ucecEst$ImmuneScore[rownames(ucecEst) %in% rownames(suhasEndo)],
        ucecEst$ImmuneScore[rownames(ucecEst) %in% rownames(suhasSerous)],
        ovEst$ImmuneScore, names = c("Endometroid EC", "Serous-like EC", "Serous OvCa"),
        ylab = "Immune score", las=2)

###
ucecCIBER = read.csv("./ucecCIBEr.csv", header = T, row.names = 1)
ovCIBER = read.csv("./ovCibersort.csv", header = T, row.names = 1)

boxplot(ucecCIBER$T.cells.CD8[rownames(ucecCIBER) %in% rownames(suhasEndo)],
        ucecCIBER$T.cells.CD8[rownames(ucecCIBER) %in% rownames(suhasSerous)],
        ovCIBER$T.cells.CD8, names = c("Endometroid", "Serous", "Ovarian"), ylab = "Immune score")

boxplot(ucecCIBER$Macrophages.M2[rownames(ucecCIBER) %in% rownames(suhasEndo)],
        ucecCIBER$Macrophages.M2[rownames(ucecCIBER) %in% rownames(suhasSerous)],
        ovCIBER$Macrophages.M2, names = c("Endometroid", "Serous", "Ovarian"), ylab = "Immune score")

boxplot(ucecCIBER$T.cells.regulatory..Tregs.[rownames(ucecCIBER) %in% rownames(suhasEndo)],
        ucecCIBER$T.cells.regulatory..Tregs.[rownames(ucecCIBER) %in% rownames(suhasSerous)],
        ovCIBER$T.cells.regulatory..Tregs., names = c("Endometroid", "Serous", "Ovarian"), ylab = "Immune score")

cibersortTotal = rbind(ucecCIBER[rownames(ucecCIBER) %in% rownames(suhasEndo),],
                       ucecCIBER[rownames(ucecCIBER) %in% rownames(suhasSerous),],
                       ovCIBER)
cibersortTotal$RMSE = NULL
cibersortTotal$P.value = NULL
cibersortTotal$Pearson.Correlation = NULL
cibersortTotal$Macrophages = cibersortTotal$Macrophages.M0+cibersortTotal$Macrophages.M1+cibersortTotal$Macrophages.M2
cibersortTotal$group = rep(c("Endometroid","Serous","Ovarian"), c(sum(rownames(ucecCIBER) %in% rownames(suhasEndo)), sum(rownames(ucecCIBER) %in% rownames(suhasSerous)), nrow(ovCIBER)))
cibersortTotal4plot = cibersortTotal %>% gather(key = "cellType",value = "abundance", 1:23)
cibersortTotal4plot$group = factor(cibersortTotal4plot$group, levels = c("Endometroid","Serous","Ovarian"))
ggplot(cibersortTotal4plot, aes(x=cellType, y=abundance, colour=group))+
  geom_boxplot()+theme(axis.text.x = element_text(angle = -90))

##
cibersortTotal$NK = cibersortTotal$NK.cells.activated+cibersortTotal$NK.cells.resting
cibersortTotal_sub = cibersortTotal[,c("T.cells.CD8","NK","T.cells.regulatory..Tregs.", "Macrophages","Monocytes","group")]
colnames(cibersortTotal_sub) = c("CD8_T_Cell","NK_Cell","Treg_Cell","Macrophages","Monocytes","Group")
cibersortTotal4plot_sub = cibersortTotal_sub %>% gather(key = "Cell_Type",value = "Abundance", 1:5)
cibersortTotal4plot_sub$Group = factor(cibersortTotal4plot_sub$Group, levels = c("Endometroid","Serous","Ovarian"))
ggplot(cibersortTotal4plot_sub, aes(x=Cell_Type, y=Abundance, colour=Group))+
  geom_boxplot()+theme(axis.text.x = element_text(angle = -90))

##
ovImmuneSurvCox=list()
for(cell in colnames(cibersortTotal_sub)[1:5]){
  ovImmuneSurvCox[[cell]] = coxph(survCombinedList$OV[ovsurvsamples,"surv"] ~ as.numeric(cibersortTotal_sub[ovsurvsamples,cell]))
}
endoImmuneSurvCox=list()
for(cell in colnames(cibersortTotal_sub)[1:5]){
  endoImmuneSurvCox[[cell]] = coxph(survCombinedList$Endo[endosurvsamples,"surv"] ~ as.numeric(cibersortTotal_sub[endosurvsamples,cell]))
}
serousImmuneSurvCox=list()
for(cell in colnames(cibersortTotal_sub)[1:5]){
  serousImmuneSurvCox[[cell]] = coxph(survCombinedList$Serous[seroussurvsamples,"surv"] ~ as.numeric(cibersortTotal_sub[seroussurvsamples,cell]))
}
survByMedian(scoreVec = setNames(as.numeric(cibersortTotal_sub[seroussurvsamples,"CD8_T_Cell"]), seroussurvsamples),
             surv = data.frame(time = suhasSerous[seroussurvsamples,"overall_survival"], 
                               event = suhasSerous[seroussurvsamples,"status"], row.names = seroussurvsamples),
             plot = T)
survByMedian(scoreVec = setNames(as.numeric(cibersortTotal_sub[endosurvsamples,"CD8_T_Cell"]), endosurvsamples),
             surv = data.frame(time = suhasEndo[endosurvsamples,"overall_survival"], 
                               event = suhasEndo[endosurvsamples,"status"], row.names = endosurvsamples),
             plot = T)
survByMedian(scoreVec = setNames(as.numeric(cibersortTotal_sub[ovsurvsamples,"CD8_T_Cell"]), ovsurvsamples),
             surv = data.frame(time = suhasOV[ovsurvsamples,"overall_survival"], 
                               event = suhasOV[ovsurvsamples,"status"], row.names = ovsurvsamples),
             plot = T)


###

mcpCounterScore = MCPcounter::MCPcounter.estimate(expression = rnaTF, featuresType = "HUGO_symbols")
mcpCounterScore = as.data.frame(t(mcpCounterScore))
cor(mcpCounterScore$`CD8 T cells`, 
               ciberAb[colnames(mcpCounterScore), `T cells CD8` ,on = "Input Sample"])
# cor(mcpCounterScore$`CD8 T cells`, 
#     ciber[rownames(mcpCounterScore), `T cells CD8` ,on = "Input Sample"])
epires = EPIC::EPIC(bulk = rnaTF)
epires = epires$cellFractions

cor.test(mcpCounterScore$`CD8 T cells`, epires[rownames(mcpCounterScore),"CD8_Tcells"])
#cor.test(mcpCounterScore$`NK cells`, epires[rownames(mcpCounterScore),"NKcells"])
cor.test(mcpCounterScore$`NK cells`, ciberAb$`NK cells resting`)
meta_table_web$MCPCounter_CD8 = mcpCounterScore[meta_table_web$CASE_ID,"CD8 T cells"]
meta_table_web$MCPCounter_NK = mcpCounterScore[meta_table_web$CASE_ID,"NK cells"]
meta_table_web$CD8_RNA = as.numeric(rnaTF["CD8A",meta_table_web$CASE_ID])
meta_table_web$CD8_Protein = as.matrix(batch13Ratio_T)["CD8A",][meta_table_web$CASE_ID]

meta_table_web$GZMA_RNA = as.numeric(rnaTF["GZMA",meta_table_web$CASE_ID])
meta_table_web$GZMA_Protein = as.matrix(batch13Ratio_T)["GZMA",][meta_table_web$CASE_ID]
meta_table_web$PRF1_RNA = as.numeric(rnaTF["PRF1",meta_table_web$CASE_ID])
meta_table_web$PRF1_Protein = as.matrix(batch13Ratio_T)["PRF1",][meta_table_web$CASE_ID]
meta_table_web$JAK2_RNA = as.numeric(rnaTF["JAK2",meta_table_web$CASE_ID])
meta_table_web$JAK2_Protein = as.matrix(batch13Ratio_T)["JAK2",][meta_table_web$CASE_ID]

meta_table_web$JAK1_RNA = as.numeric(rnaTF["JAK1",meta_table_web$CASE_ID])
meta_table_web$JAK1_Protein = as.matrix(batch13Ratio_T)["JAK1",][meta_table_web$CASE_ID]


####
cor.test(meta_table_web$MCPCounter_CD8, meta_table_web$PDL1_Protein)
cor.test(meta_table_web$MCPCounter_CD8, meta_table_web$PDL1_RNA)
cor.test(meta_table_web$MCPCounter_CD8, meta_table_web$PDL1_RNA)
cor.test(meta_table_web$CD8_RNA, meta_table_web$PDL1_RNA)
cor.test(meta_table_web$CD8_Protein, meta_table_web$PDL1_RNA)

appGeneList = read.table("~/Downloads/kegg_app.txt", stringsAsFactors = F)
appGeneList = list(app = appGeneList$V1)

###
appScoreRNA = GSVA::gsva(expr = as.matrix(rnaTF), gset.idx.list = appGeneList, method = "ssgsea")
appScoreProtein = GSVA::gsva(expr = as.matrix(batch13Ratio_T), 
                             gset.idx.list = appGeneList, method = "ssgsea")

cor(appScoreRNA[1, colnames(appScoreProtein)],appScoreProtein[1,])

meta_table_web$APP_score_RNA = appScoreRNA[1, meta_table_web$CASE_ID]
meta_table_web$APP_score_Protein = appScoreProtein[1, ][meta_table_web$CASE_ID]
meta_table_web$APP_score_RNA = appScoreRNA[1, meta_table_web$CASE_ID]


meta_table_web$chrIdx = colSums(chrIndex)[meta_table_web$CASE_ID]
meta_table_web$chr3pLoss = arm_loss["chr3p",meta_table_web$CASE_ID]
meta_table_web$chr3qGain = arm_gain["chr3q",meta_table_web$CASE_ID]
meta_table_web$chr9pLoss = arm_loss["chr9p",meta_table_web$CASE_ID]
meta_table_web$chr11pGain = arm_gain["chr11q",meta_table_web$CASE_ID]
quickImmuneMatrix = meta_table_web[order(CD8_RNA),
                                   .(chrIdx,chr3pLoss,chr3qGain, chr9pLoss, chr11pGain,`A.Mutation Count`,JAK1_RNA, JAK1_Protein, JAK2_RNA, JAK2_Protein,JAK.STAT, PDL1_RNA, PDL1_Protein, APP_score_RNA, 
                                     APP_score_Protein, CD8_RNA, CD8_Protein,
                                     GZMA_RNA,GZMA_Protein,PRF1_RNA,PRF1_Protein, 
                                     MCPCounter_CD8, MCPCounter_NK)]

quickImmuneMatrix = t(quickImmuneMatrix)                         
quickImmuneMatrix = t(scale(t(quickImmuneMatrix)))
colnames(quickImmuneMatrix) = meta_table_web[order(CD8_RNA),CASE_ID]

quickImmuneAnno = meta_table_web[,.(TUMOR_SITE_CURATED, TUMOR_STAGE_CURATED, 
                                    tanscriptomic_subtype_new, Smoking.Influence)]
quickImmuneAnno = data.frame(quickImmuneAnno)
rownames(quickImmuneAnno) = meta_table_web$CASE_ID

Heatmap(quickImmuneMatrix, cluster_rows = F, cluster_columns = F, 
        name = "scaled_value", 
        top_annotation = HeatmapAnnotation(df =quickImmuneAnno[colnames(quickImmuneMatrix),],
                                           col = list(TUMOR_SITE_CURATED=setNames(brewer.pal(5,"Accent"), c("Oral cavity", "Oropharynx", "Larynx","Lip","Hypopharynx")),
                                                      TUMOR_STAGE_CURATED = colorDiscreetCodeFunc(vecSeq = 4,start = "yellow",col = "purple",vecName = paste0("Stage ",c("I","II","III","IV"))),
                                                      tanscriptomic_subtype_new = setNames(brewer.pal(4,"Set3"), c("Basal","Classical","Atypical","Mesenchymal")),
                                                      Smoking.Influence = colorDiscreetCodeFunc(vecSeq = 3, vecName = c("cancer without smoking evidence","cancer with minor/limited smoking evidence","cancer with smoking evidence"))
                                                      )),
        show_column_names = F)

ctGenes = c("KIF2C","CEP55", "NUF2", "TTK", "SPAG9", "MIA2", "MAGEA1", "KDM5B", "CNOT9", "IGF2BP3")
ctAntigene = batch13Ratio_T[ctGenes,]      
ctAntigene = cbind(ctAntigene, matrix(NAs(length(need2fill)*length(ctGenes)), nrow = length(ctGenes)))
colnames(ctAntigene) = c(colnames(batch13Ratio_T), need2fill)
ctAntigene = rbind(ctAntigene, rnaTF[ctGenes, colnames(ctAntigene)])
ctAntigene = t(scale(t(ctAntigene)))

rownames(ctAntigene) = paste0(rownames(ctAntigene), "_Protein")
rownames(ctAntigene)[10:20] = gsub(rownames(ctAntigene)[10:20], pattern = "1_Protein", replacement = "_RNA")
ctAntigene = t(scale(t(ctAntigene)))
Heatmap(ctAntigene[,colnames(quickImmuneMatrix)], cluster_rows = F, cluster_columns = F, name = "scaled_value_CT_antigen")


#meta_table_web[,1:ncol(meta_table_web):=lapply(.SD, function(x){x[x=="unknown"]=NA; x}),]
#meta_table_web[,1:ncol(meta_table_web):=lapply(.SD, function(x){x[x=="NA"]=NA; x}),]
                      