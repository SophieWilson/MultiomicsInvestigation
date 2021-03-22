library(mixOmics)
setwd('C:/Users/Mischa/Documents/Uni Masters/Module 5- complex')

# loading the data
rna_df <- read.csv('C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/rna_vst_counts.csv', header = TRUE, row.names = 1)
sample_sheet <- read.csv('C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/sample_sheet.csv')

# formatting the data
rna_class <- as.factor(sample_sheet$REF)
summary(rna_class)
rna_dft <- t(rna_df)
dim(rnadf)
length(rna_class)


####### Tuning the hyperparameters of SPLSDA 
# looking for ncomp to include
hyperparameter_tune <- function(plsda_result, class_list, dataframe){
  set.seed(30) # for reproducbility in this vignette, otherwise increase nrepeat
  MyPerf.plsda <- perf(plsda_result, validation = "Mfold", folds = 3, 
                       progressBar = FALSE, nrepeat = 50) # we suggest nrepeat = 50
  list.keepX <- c(5:10,  seq(20, 100, 10)) # then seeing how much to keepX for each
  tune.splsda.srbct <- tune.splsda(dataframe, class_list, ncomp = 5, validation = 'Mfold', folds = 3, dist = 'max.dist', progressBar = FALSE, measure = "BER", test.keepX = list.keepX, nrepeat = 50)   # we suggest nrepeat = 50
  error <- tune.splsda.srbct$error.rate # error rate
  ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests on the error rate
  select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
  MyResult.splsda.final <- splsda(dataframe, class_list, ncomp = ncomp, keepX = select.keepX)
  
  tuning_results <- list(perf_plsda = MyPerf.plsda, keepXtested = list.keepX, error = error, ncomp = ncomp, bestkeepX = select.keepX, final = MyResult.splsda.final, tune.spla.srbct = tune.splsda.srbct)
  return(tuning_results)
}

MyResult.plsda2 <- plsda(rna_dft,rna_class, ncomp=10)
auc.untunedRNA <- auroc(MyResult.plsda2)
RNA_hyper <- hyperparameter_tune(MyResult.plsda2, rna_class, rna_dft)

#### plotting tuning results


plot(RNA_hyper$perf_plsda, col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")
plot(RNA_hyper$tune.spla.srbct)
plotIndiv(RNA_hyper$final, ind.names = FALSE, legend=TRUE, ellipse = TRUE, title="sPLS-DA - final result")

##### plotting results

results.splsda <- RNA_hyper$final

# plotting
plotIndiv(results.splsda)
# seperation between some groups but not control
plotVar(results.splsda)
# no differentiation between groups on either component
## making some better graphs
#USE THIS ONE
plotIndiv(results.splsda, ind.names = FALSE, legend=TRUE, ellipse = TRUE, star = TRUE, title = 'sPLS-DA on RNA', X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')
plotVar(results.splsda, var.names=FALSE)
# This still shows nothing, everything correlated on both components. 
#If you set the cutoff to 0.7 you can see some important variables. 
plotVar(results.splsda, cutoff = 0.7)
# plotting an auc to see how good it was at predicting the groups
auc.plsda <- auroc(results.splsda)
#Variable selection
x <- selectVar(results.splsda, comp=1)$value
x$gene <- rownames(x)
# getting the top 20 genes that contribute to component 1
top_20_genesSPLSDAc1 <- head(x[order(-x$value.var),], n=20)
# doing the same for component 2
x2 <- selectVar(results.splsda, comp=2)$value
x2$gene <- rownames(x2)
top_20_genesSPLSDAc2 <- head(x2[order(-x2$value.var),], n=20)
# plotting all
plotLoadings(results.splsda, contrib = 'max', method = 'mean')
# plotting only top 20
plotLoadings(results.splsda, contrib = 'max', method = 'mean', ndisplay = 20)
#comp 2
plotLoadings(results.splsda, contrib = 'max', method = 'mean', comp = 2)




##############################################################################

# Doing the same for metabolomics data. 
# reading it in
metab_df <- read.csv('polar_pos_pqn_imputed_glog.csv', row.names = 1)
# transposing
metab_df <- t(metab_df)
# removing the mz row
metab_df <- metab_df[-1,]
# load sample sheet
sample_sheet <- read.csv('C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/sample_sheet.csv')
# first ordering the sheet by sample name (just in case it isnt)
metab_sheet <- sample_sheet[order(sample_sheet$SampleID),]
# renumbering the rows 
rownames(metab_sheet) <- 1:nrow(metab_sheet)
# removing the sample missing from the metabolomics df
metab_sheet <- metab_sheet[-16,]
metab_class <- as.factor(metab_sheet$REF)


# doing a sparse plsda and plotting the results
metabresults.splsda2 <- splsda(metab_df, metab_class, ncomp = 10)
metauntuned <- auroc(metabresults.splsda2)
metab_hyper <- hyperparameter_tune(metabresults.splsda, metab_class, metab_df) 
metabresults.splsda <- metab_hyper$final

plotIndiv(metab_hyper$final)
# seperation between some groups but not control
plotVar(metab_hyper$final)
# no differentiation between groups on either component

## making some better graphs
#USE THIS ONE
plotIndiv(metabresults.splsda, ind.names = FALSE, legend=TRUE, ellipse = TRUE, star = TRUE, title = 'sPLS-DA on Metabolomics', X.label = 'PLS-DA 1', Y.label = 'PLS-DA 2')

plotVar(metabresults.splsda, var.names=FALSE)
# This still shows nothing, everything correlated on both components. 
#If you set the cutoff to 0.7 you can see some important variables. 
plotVar(metabresults.splsda, cutoff = 0.6)
# should figure out how to get a list of these variables
# plotting an auc to see how good it was at predicting the groups
auc.plsdameta <- auroc(metabresults.splsda)


#Variable selection
xmeta <- selectVar(metabresults.splsda, comp=1)$value
xmeta$mz <- rownames(xmeta)
# getting the top 20 genes that contribute to component 1
top_20_mzSPLSDAc1 <- head(xmeta[order(-xmeta$value.var),], n=20)

# doing the same for component 2
xmeta2 <- selectVar(metabresults.splsda, comp=2)$value
xmeta2$mz <- rownames(xmeta2)
top_20_mzSPLSDAc2 <- head(xmeta2[order(-xmeta2$value.var),], n=20)

# trying to plot these
# comp1
# plotting all
plotLoadings(metabresults.splsda, contrib = 'max', method = 'mean')
# plotting only top 20
plotLoadings(metabresults.splsda, contrib = 'max', method = 'mean', ndisplay = 20)
#comp 2
plotLoadings(metabresults.splsda, contrib = 'max', method = 'mean', ndisplay = 20, comp = 2)

# Writing out the results (the tuned results)
# RNA top genes
write.csv(top_20_genesSPLSDAc1, 'C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/PLSDA output/top_20_genesPLSDAc1.csv')
write.csv(top_20_genesSPLSDAc2, 'C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/PLSDA output/top_20_genesPLSDAc2.csv')
# metabolome to mz peaks
write.csv(top_20_mzSPLSDAc1, 'C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/PLSDA output/top_20_mzSPLSDAc1.csv')
write.csv(top_20_mzSPLSDAc2, 'C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/PLSDA output/top_20_mzSPLSDAc2.csv')



