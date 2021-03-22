library(DESeq2)
library(edgeR)
library(DEFormats)

## RNA seq analysis
setwd('C:/Users/Mischa/Documents/Uni Masters/Module 5- complex')
rna_norm_df <- read.csv('rna_raw_counts.csv', header = TRUE, row.names = 1)
sample_sheet <- read.csv('sample_sheet.csv')


sample_sheet2 <- sample_sheet[-c(1:6),]

rna_norm_df2 <- rna_norm_df[,colnames(rna_norm_df) %in% sample_sheet2$SampleID ]

dge <- DGEList(round(rna_norm_df2), group = sample_sheet2$REF)
dds <- as.DESeqDataSet(dge)
sdd <- DESeq(dds)
res <- results(sdd)
resOrdered <- res[order(res$pvalue),]
summary(res)
res05 <- results(sdd)
summary(res05)
plotMA(res)


# shrinking log fold change for visualisation and ranking
resultsNames(sdd)

resLFC2 <- lfcShrink(sdd, coef = 'group_1x_vs_10x', type='apeglm')
plotMA(resLFC2, ylim = c(-1,2))
summary(resLFC2)
sum(res$padj < 0.05, na.rm = TRUE)

df <- data.frame(rownames(res), res$padj)
df <- na.omit(df)
sig = ''
for (i in 1:length(df$res.padj)){
  if (df$res.padj[i] < 0.05){
    append(sig,(df$rownames.res[i]))
  }
}

imp_MA <- ''
any(is.na(resLFC$pvalue))
length(resLFC$pvalue)
no <- na.omit(resLFC2$pvalue)
length(no)
for (i in 1:length(no)){
  if (no[i] < 0.1){
    imp_MA[i] = no[i]
  }
}
any(is.na(res$pvalue)) 
indx <- apply(res$pvalue, 1, function(x) any(is.na(x)))

### Water chemical analysis

# loading the water chemicals
water_chemicals <- read.table('water_chemicals.tsv',sep = '\t', header = TRUE, row.names = 1)
# removing the chemicals with all 0 values
water2 <- water_chemicals[rowSums(water_chemicals == 0) <= 11, ]
write.csv(water2, file = 'waterchemicalsnozero.csv')
# making a heatmap, all the same
cordf <- cor(water2[,-1])
library(reshape2)
melted_cordf <- melt(cordf)
library(ggplot2)
p <- ggplot(data = melted_cordf, aes(x=X1, y=X2, fill=value)) + 
  geom_tile() 
p + labs(color = "Cylinders\n", x = 'Var1', y = 'Var2') + ggtitle('Stream Location Composition Correlation')


# trying some line charts who knows
library(reshape)
water2$chemical <- rownames(water2)
Molten <- melt(water2, id.vars = 'chemical')
Molten <- Molten[-c(1:16),]
ggplot(Molten, aes(x = variable, y = value, group = chemical)) + geom_line(aes(color = chemical))
# this plot showed no overall trends

# heatmap
ggplot(Molten, aes(x = variable, y = value, fill = value)) + geom_tile()

## i guess ill try an anova? 
library(multcomp)
Molten$value <- as.numeric(Molten$value)
res.aov <- aov(value ~ variable, data = Molten)
st <- summary(glht(res.aov, linfct=mcp(variable='Tukey')$'PR(>|t|)'))
save(st, file = 'anovawaterchemicaldiff')
## No anova was signficiant, no significant difference between the groups (locations)

# tried it for individual chemicals, it didnt work. will assume no relationship
test <- melt(water2, id.vars = 'chemical')
Moltenchem <- Molten[order(Molten$chemical),]
Moltenchem$value <-Moltenchem$value + 0.001
er <- Moltenchem[c(1:11),]
res2.aov <- aov(value ~ variable, data = er)
summary(glht(res2.aov, linfct=mcp(variable='Tukey')))


## PCA on RNA data?
rna_norm_tsv <- read.csv('rna_norm_counts.csv', header = TRUE, row.names = 1)
# # transforming it so that samples are on the row and chemicals are columns
pca_x <- t(rna_norm_tsv)
# # doing the pca
pca_res1 <- prcomp(pca_x, center = TRUE, scale = TRUE)
summary(pca_res1)
# making the dataframe
pc_df1 <- data.frame(pca_res1$x, Condition = sample_sheet$REF, Location = sample_sheet$Site)
# plotting it

library(RColorBrewer)
nb.cols <- length(unique(pc_df1$Location))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
ggplot(pc_df1, aes(PC1, PC2, shape = Condition)) + 
  geom_point(aes(color = Location), size = 5) + scale_color_manual(values = mycolors)


# PCA for metabolomics
metabo_pos <- read.csv('polar_pos_pqn_imputed_glog.csv', header = TRUE, row.names = 1)
metabo_pos <- metabo_pos[,-1]
# # transforming it so that samples are on the row and chemicals are columns
setdiff(sample_sheet$SampleID, colnames(metabo_pos))
# removing the pesky sample
sample_meta <- sample_sheet[-82,]
pca_meta <- t(metabo_pos)
# # doing the pca
pca_res2 <- prcomp(pca_meta, center = TRUE, scale = TRUE)
summary(pca_res2)
# making the dataframe
pc_df2 <- data.frame(pca_res2$x, Condition = sample_meta$REF, Location = sample_meta$Site)
# plotting it
ggplot(pc_df2, aes(x = PC1, y = PC2, shape = Condition)) + geom_point(aes(color = Location), size = 5) + scale_color_manual(values = mycolors)


hmdbres <- read.table(file = 'Annotation/polar_pos_peaklist_hmdb_results.tsv', sep = '\t', header = TRUE)
keggres <- read.table(file = 'Annotation/polar_pos_peaklist_kegg_results.tsv', sep = '\t', header = TRUE)

write.csv(hmdbres, 'C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/Annotation/hmdbres.csv')

write.csv(keggres, 'C:/Users/Mischa/Documents/Uni Masters/Module 5- complex/Annotation/keggeres.csv')
