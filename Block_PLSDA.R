library(mixOmics)
setwd('C:/Users/Mischa/Documents/Uni Masters/Module 5- complex')

# loading the data
rna_df <- read.csv('C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/rna_vst_counts.csv', header = TRUE, row.names = 1)
sample_sheet <- read.csv('C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/sample_sheet.csv')
metab_df <- read.csv('polar_pos_pqn_imputed_glog.csv', row.names = 1)


# formatting the data
rna_dft <- t(rna_df)
rna_dft <- rna_dft[-82,] # removing the pesky sample
rna_dft <- rna_dft[order(rownames(rna_dft)),]
metab_df <- metab_df[order(rownames(metab_df)),]
# Metabolomics data. 
metab_df <- t(metab_df)
metab_df <- metab_df[-1,] # removing the mz row
metab_sheet <- sample_sheet[order(sample_sheet$SampleID),]
metab_sheet <- sample_sheet[-82,] # removing the sample missing from the df

metab_class <- as.factor(metab_sheet$REF)
# setting variables
X <- list(Rna = rna_dft, Metabolomics = metab_df)
Y <- metab_class
MyResult.diablo <- block.splsda(X, Y, ncomp = 10)
auroc(MyResult.diablo)
plotLoadings(MyResult.diablo, comp = 1, contrib = "max", ndisplay = 150)
plotIndiv(sgccda.res)


# Tuning hyperparameters
set.seed(123)
MyPerf.diablo <- perf(MyResult.diablo, validation = 'Mfold', folds = 6, nrepeat = 10, dist = 'centroids.dist')
perf.diablo$choice.ncomp$WeightedVote
plot(perf.diablo)

design = matrix(0.1, ncol = length(X), nrow = length(X), 
                dimnames = list(names(X), names(X)))
diag(design) = 0

design 

ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] # keeping the optimal ncomp
test.keepX = list (Rna = c(5:10,  seq(20, 10)),
                   Metabolomics = c(5:10,  seq(20, 10)))

tune.TCGA = tune.block.splsda(X = X, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design, validation = 'Mfold', folds = 10, nrepeat = 1, cpus = 2, dist = "centroids.dist")
list.keepXmy <- list(Rna = c(2,10), Metabolomics = c(5,10))
list.keepX = tune.TCGA$choice.keepX # tuning keepX
# the final model
sgccda.res = block.splsda(X = X, Y = Y, ncomp = 2, 
                          keepX = list.keepXmy, design = design)


perfblock <- auroc(sgccda.res)
saveRDS(sgccda.res, file="perfblockparams.RData")
test <- readRDS("C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/perfblockparams.RData")



plotIndiv(MyResult.diablo, ind.names = FALSE, legend=TRUE) ## sample plot
plotVar(sgccda.res) ## variable plot
plotIndiv(sgccda.res, ind.names = FALSE, legend=TRUE,ellipse = TRUE, star = TRUE, cex=c(1,2,3), title = 'Block SPLSDA', X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')
plotIndiv(MyResult.diablo, ind.names = FALSE, legend=TRUE, ellipse = TRUE, star = TRUE, title = 'sPLS-DA on RNA', X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')
plotVar(sgccda.res, var.names = c(FALSE, FALSE), legend=TRUE, pch=c(5,5))
plotDiablo(sgccda.res, ncomp = 1)

circos<- circosPlot(sgccda.res, cutoff = 0.7,showIntraLinks = TRUE, line = TRUE, size.variables = 0.61, size.labels = 1.4, var.adj = -0.5)
SelectedVariablescircos <- selectVar(sgccda.res, comp =2)
SelectedVariablescircos$Rna$name
SelectedVariablescircos$Metabolomics$name


cimDiablo(sgccda.res, margin=c(8,20))
cimDiablo(sgccda.res, color.blocks = c('darkorchid', 'lightgreen'), comp = 2, margin=c(8,20), legend.position = "right")
plotLoadings(sgccda.res, comp = 2, contrib = "max")
#network(MyResult.diablo, blocks = c(1,2,3),color.node = c('darkorchid', 'brown1', 'lightgreen'), cutoff = 0.6, save = 'jpeg', name.save = 'DIABLOnetwork')

auc.blcokplsda <- auroc(MyResult.diablo)





#Variable selection RNA
#c1
varselectRNA <- selectVar(MyResult.diablo, comp=1)$rna$value
varselectRNA$gene <- rownames(varselectRNA)
# getting the top 20 genes that contribute to component 1
RNA_blockplsda_c1_20 <- head(varselectRNA[order(-varselectRNA$value.var),], n=20)
#c2
varselectRNA2 <- selectVar(MyResult.diablo, comp=2)$rna$value
varselectRNA2$gene <- rownames(varselectRNA2)
RNA_blockplsda_c2_20 <- head(varselectRNA2[order(-varselectRNA2$value.var),], n=20)

# Variable selection Metabolomics
#c1
varselectMETA <- selectVar(MyResult.diablo, comp=1)$metabolomics$value
varselectMETA$mz <- rownames(varselectMETA)
# getting the top 20 genes that contribute to component 1
META_blockplsda_c1_20 <- head(varselectMETA[order(-varselectMETA$value.var),], n=20)
#c2
varselectMETA2 <- selectVar(MyResult.diablo, comp=2)$metabolomics$value
varselectMETA2$mz <- rownames(varselectMETA2)
# getting the top 20 genes that contribute to component 1
META_blockplsda_c2_20 <- head(varselectMETA2[order(-varselectMETA$value.var),], n=20)

# plotting all
#rna
png(file="testplot.png")
plotLoadings(MyResult.diablo, contrib = 'max', method = 'mean')
dev.off()

# plotting only top 20
plotLoadings(MyResult.diablo, contrib = 'max', method = 'mean', ndisplay = 20)
#comp 2
plotLoadings(MyResult.diablo, contrib = 'max', method = 'mean', ndisplay = 20, comp = 2)

# writing out the files
# RNA 
write.csv(RNA_blockplsda_c1_20, 'C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/Block PLSDA/RNA_blockplsda_c1_20.csv')
write.csv(RNA_blockplsda_c2_20, 'C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/Block PLSDA/RNA_blockplsda_c2_20.csv')
# metabolome to mz peaks
write.csv(META_blockplsda_c1_20, 'C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/Block PLSDA/META_blockplsda_c1_20.csv')
write.csv(META_blockplsda_c2_20, 'C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/Block PLSDA/META_blockplsda_c2_20.csv')

