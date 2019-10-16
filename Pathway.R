#This file contain all the scripts used for pathway-related analysis

##progeny pathways #not run
library(progeny)
all(rownames(ucecRNA) == rownames(ovRNA))
ucecPath = progeny(expr = as.matrix(cbind(ucecRNA,ovRNA)), scale = TRUE)
hallMark = qusage::read.gmt("~/mkProj/h.all.v6.2.symbols.gmt")
library(GSVA)
hallMarkScore = gsva(expr = as.matrix(cbind(ucecRNA,ovRNA)), gset.idx.list = hallMark)
rnaAnno = data.frame(ID = c(intersect(colnames(ucecRNA), rownames(suhasEndo)),
                            intersect(colnames(ucecRNA), rownames(suhasSerous)),
                            colnames(ovRNA)),
                     type = rep(c("Endo","Serous","OV"),c(395, 110, 303)),
                     stringsAsFactors = F)
rownames(rnaAnno) = rnaAnno$ID
hallMarkScore = hallMarkScore[,rnaAnno$ID]
autoplot(prcomp(t(hallMarkScore)), data = rnaAnno, colour = "type", frame = TRUE, frame.type = 'norm')
hallMarkScore = t(hallMarkScore)
hallMarkScore = as.data.frame(hallMarkScore)
hallMarkScore$type = rnaAnno$type
hallMarkScore4plot = hallMarkScore %>% gather(value = "score", key  = "pathway", 1:50)
hallMarkScore4plot$type = factor(hallMarkScore4plot$type, levels = c("Endo","Serous","OV"))
or
ggplot(hallMarkScore4plot, aes(x = pathway, y = score, color=type))+geom_boxplot()

##only show anova siganifcant results
anovaList = list()
for(path in colnames(hallMarkScore)[1:50]){
  anovaList[[path]] = aov(hallMarkScore[,path] ~ hallMarkScore$type)
}
anovaP = sapply(anovaList, function(x) summary(x)[[1]][["Pr(>F)"]])
anovaP = anovaP[1,]
anovaPadj = p.adjust(anovaP)
table(anovaPadj<0.01)
hallMarkScoreMedian = aggregate(hallMarkScore[,1:50], by=list(hallMarkScore$type), FUN=median)
rownames(hallMarkScoreMedian) = hallMarkScoreMedian$Group.1
hallMarkScoreMedian$Group.1 = NULL
hallMarkScoreMedian = t(hallMarkScoreMedian)
hallMarkScoreMedian = hallMarkScoreMedian[,c(1,3,2)]
sigAnno = rowAnnotation(df = data.frame(row.names = names(anovaPadj), significance = anovaPadj<0.01))
Ht = ComplexHeatmap::Heatmap(hallMarkScoreMedian, clustering_distance_rows = "spearman",cluster_columns = F)
#   clustering_distance_columns = "euclidean")
draw(sigAnno+Ht)

hallMarkScoreMedian_sig = hallMarkScoreMedian[names(anovaPadj<0.01),]
rownames(hallMarkScoreMedian_sig) = gsub(rownames(hallMarkScoreMedian_sig), pattern = "HALLMARK_", replacement = "")
ComplexHeatmap::Heatmap(hallMarkScoreMedian_sig, clustering_distance_rows = "spearman", 
                        name = "Pathway\nActivity", row_names_gp = grid::gpar(fontsize = 5), cluster_columns = F)
#cluster_columns = T, clustering_distance_columns = "euclidean")
