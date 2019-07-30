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