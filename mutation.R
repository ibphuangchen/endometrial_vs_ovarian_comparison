#this file includes the scripts for somatic mutation analysis
ovMut = read.table("./OV/Human__TCGA_OV__WUSM__Mutation__GAIIx__01_28_2016__BI__Gene__Firehose_MutSig2CV.cbt",
                   sep = "\t",header = T, stringsAsFactors = F, row.names = 1)
dim(ovMut)
endoMut = read.table("./OV/Human__TCGA_UCEC__WUSM__Mutation__GAIIx__01_28_2016__BI__Gene__Firehose_MutSig2CV.cbt",
                     sep = "\t",header = T, stringsAsFactors = F, row.names = 1)
dim(endoMut)
rm(endoMut)

library(TCGAbiolinks)
ovMaf <- GDCquery_Maf("OV", pipelines = "mutect")
length(unique(ovMaf$Tumor_Sample_Barcode)) #436

ucecMaf = readRDS("~/ucecTCGA/ucecMaf.rds")
length(unique(ucecMaf$Tumor_Sample_Barcode)) #530
library(dplyr)
ucecMaf = ucecMaf %>% filter(grepl(Tumor_Sample_Barcode, pattern = "01A"))
ucecMaf$Tumor_Sample_Barcode = gsub(ucecMaf$Tumor_Sample_Barcode, pattern = "01A.*", replacement = "")
ovMaf = ovMaf %>% filter(grepl(Tumor_Sample_Barcode, pattern = "01A"))
ovMaf$Tumor_Sample_Barcode = gsub(ovMaf$Tumor_Sample_Barcode, pattern = "01A.*", replacement = "")
endoMaf = ucecMaf %>% filter(Tumor_Sample_Barcode %in% names(ucecHis)[ucecHis == "Endometrioid"])
length(unique(endoMaf$Tumor_Sample_Barcode)) #391
SerousMaf = ucecMaf %>% filter(Tumor_Sample_Barcode %in% names(ucecHis)[ucecHis == "Serous"])
length(unique(SerousMaf$Tumor_Sample_Barcode)) #110
SerousMutStat = SerousMaf %>% group_by(Hugo_Symbol) %>% distinct(Tumor_Sample_Barcode) %>% count(Hugo_Symbol,name = "count") %>% arrange(desc(count))
SerousMutStat$proportion = SerousMutStat$count/110
EndoMutStat = endoMaf %>% group_by(Hugo_Symbol) %>% distinct(Tumor_Sample_Barcode) %>% count(Hugo_Symbol,name = "count") %>% arrange(desc(count))
EndoMutStat$proportion = EndoMutStat$count/391
ovMutStat = ovMaf %>% group_by(Hugo_Symbol) %>% distinct(Tumor_Sample_Barcode) %>% count(Hugo_Symbol,name = "count") %>% arrange(desc(count))
ovMutStat$proportion = ovMutStat$count/436
mutPlot = rbind(SerousMutStat[1:20,],EndoMutStat[1:20,],ovMutStat[1:20,])
mutPlot$type = rep(c("Serous","Endometroid","Ovarian"),each = 20)
mutPlot$proportion = mutPlot$count/rep(c(110,391,436), each=20)
mutPlot = split(mutPlot, f=mutPlot$type)
mutPlot = lapply(mutPlot, function(x)x[x$Hugo_Symbol!="TTN",])
plotList = list()
for(i in 1:3){
  mutPlot[[i]]$Hugo_Symbol = factor(mutPlot[[i]]$Hugo_Symbol, levels = mutPlot[[i]]$Hugo_Symbol)
  plotList[[i]]=ggplot(mutPlot[[i]], aes(x=Hugo_Symbol,y=proportion))+
    geom_bar(stat = "identity")+xlab("")+ylab("")+ylim(0,1)+
    theme(axis.text.x = element_text(angle = -90))
}
plotList
rm(ucecMaf)
rm(ovMaf)
rm(endoMaf)
rm(SerousMaf)


####calculate the mutation load
endoCount = endoMaf %>% count(Tumor_Sample_Barcode,name = "mutation_load")
SerousCount = SerousMaf %>% count(Tumor_Sample_Barcode,name = "mutation_load")
ovCount = ovMaf %>% count(Tumor_Sample_Barcode,name = "mutation_load")
combineCount = rbind(endoCount, SerousCount, ovCount)
combineCount$type = rep(c("Endometroid EC","Serous-like EC","serous OvCa"),c(nrow(endoCount), nrow(SerousCount), nrow(ovCount)))
combineCount$type = factor(combineCount$type, levels = c("Endometroid EC","Serous-like EC","serous OvCa"))
ggplot(combineCount, aes(x=type, y=mutation_load))+
  geom_boxplot(outlier.shape = NA)+ylim(0,2000)+
  ylab("Mutation load score")+xlab("")+
  theme_classic()+theme(axis.text.x = element_text(angle = -90, size = 10))
  
t.test(endoCount$mutation_load, SerousCount$mutation_load)

##comparison for specific sfes: Ciriello et al, 2015
sfe = read.csv("./sfe.csv", stringsAsFactors = F)
sfe = sfe$Alteration
sfe = sfe[!grepl(sfe, pattern = "chr")]
length(sfe)
sfe = intersect(sfe,rownames(endoGeneLevel))

EndoMutStat= as.data.frame(EndoMutStat); rownames(EndoMutStat) = EndoMutStat$Hugo_Symbol
SerousMutStat= as.data.frame(SerousMutStat); rownames(SerousMutStat) = SerousMutStat$Hugo_Symbol
ovMutStat= as.data.frame(ovMutStat); rownames(ovMutStat) = ovMutStat$Hugo_Symbol

endoGeneLevel_scna_mut = data.frame(scna = apply(endoGeneLevelB[sfe,], 1, function(x)sum(abs(x)==2)/length(x)),
                                    mutation = EndoMutStat[sfe,"proportion"])
serousGeneLevel_scna_mut = data.frame(scna = apply(serousGeneLevelB[sfe,], 1, function(x)sum(abs(x)==2)/length(x)),
                                    mutation = SerousMutStat[sfe,"proportion"])
serousGeneLevel_scna_mut$mutation[is.na(serousGeneLevel_scna_mut$mutation)] = 0
ovGeneLevel_scna_mut = data.frame(scna = apply(ovGeneLevelB[sfe,], 1, function(x)sum(abs(x)==2)/length(x)),
                                    mutation = ovMutStat[sfe,"proportion"])
ovGeneLevel_scna_mut$mutation[is.na(ovGeneLevel_scna_mut$mutation)] = 0
ggplot(endoGeneLevel_scna_mut, aes(scna,mutation))+geom_point()
ggplot(serousGeneLevel_scna_mut, aes(scna,mutation))+geom_point()
ggplot(ovGeneLevel_scna_mut, aes(scna,mutation))+geom_point()
wilcox.test(endoGeneLevel_scna_mut$scna, endoGeneLevel_scna_mut$mutation, alternative = "less")
wilcox.test(serousGeneLevel_scna_mut$scna, endoGeneLevel_scna_mut$mutation, alternative = "less")
wilcox.test(ovGeneLevel_scna_mut$scna, endoGeneLevel_scna_mut$mutation, alternative = "less")
wilcox.test(endoGeneLevel_scna_mut$scna, endoGeneLevel_scna_mut$mutation, alternative = "greater")
wilcox.test(serousGeneLevel_scna_mut$scna, endoGeneLevel_scna_mut$mutation, alternative = "greater")
wilcox.test(ovGeneLevel_scna_mut$scna, endoGeneLevel_scna_mut$mutation, alternative = "greater")
