rm(list=ls())
library("limma")
library("edgeR")

#read in raw read counts (RNA-seq data available from GEO)
x <- read.table(file="data/Koch_Counts.txt", header=T, check.names = F, row.names=1)

#sample names and information (selection regime, treatment, batch=sequencing run)
samples.Transplant.seq<- read.table(file="data/Samples_Transpl_2018.txt",header=T)

#select samples in CT and HD
CT.HD.index <- which(samples.Transplant.seq$group=="CT.CT"|
                      samples.Transplant.seq$group=="HD.CT"|
                      samples.Transplant.seq$group=="CT.HD"|
                      samples.Transplant.seq$group=="HD.HD")

counts.CT.HD <- x[,CT.HD.index]
Samples.CT.HD <- samples.Transplant.seq[CT.HD.index,]
##### !Remove extreme outlier Mx1CT-1-4-CT (probably mistake during library preparation)
Samples.CT.HD <- Samples.CT.HD[which(!Samples.CT.HD$Sample=="M1CT-1-4CT"),]
counts.CT.HD <- counts.CT.HD[, which(!colnames(counts.CT.HD)=="Mx1CT-1-4CT")]


Samples.CT.HD$Line <- paste(substr(Samples.CT.HD$Sample,1,2), Samples.CT.HD$Selection, sep="")
Samples.CT.HD$group <- factor(Samples.CT.HD$group, levels=c("CT.CT", "HD.CT", "CT.HD", "HD.HD")) 
Samples.CT.HD$Line <- factor(Samples.CT.HD$Line)
Samples.CT.HD$Selection <- factor(Samples.CT.HD$Selection)

y <- DGEList(counts=counts.CT.HD,group=Samples.CT.HD$group)
design <- model.matrix(~0 + batch+group, data=Samples.CT.HD)

#filter lowly expressed genes
keep <- rowSums(cpm(y)>1) >=2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- DGEList(y,group=Samples.CT.HD$group)
y <- calcNormFactors(y)
v <- voom(y, design)
# accounting for non-independence of samples: biological replicates (selection lines)
#and technical replicates (samples from the same selection line)
corfit <- duplicateCorrelation(v, design, block = Samples.CT.HD$Line)
v <- voom(y, design, block = Samples.CT.HD$Line, correlation =
            corfit$consensus)
corfit <- duplicateCorrelation(v, design, block = Samples.CT.HD$Line)

fit <- lmFit(v, design, block = Samples.CT.HD$Line, correlation =
               corfit$consensus)
##making contrasts to test for differences in baseline expression in CT and HD
## and differences in responses to HD treatment
baselineDiff.inCT <- makeContrasts( groupHD.CT,
                                    levels=design)
baselineDiff.inHD <- makeContrasts( groupHD.HD-groupCT.HD,
                                   levels=design)
responseCT <- makeContrasts( groupCT.HD,
                             levels=design)
responseHD <- makeContrasts( groupHD.HD-groupHD.CT,
                            levels=design)
total.change <- makeContrasts( groupHD.HD,
                               levels=design)
Diff.responses <- makeContrasts((groupHD.HD-groupHD.CT)-(groupCT.HD),
                                levels=design)

#Expression differences between Hot.Dry-line and Control-lines in control conditions 
fit2 <- contrasts.fit(fit, baselineDiff.inCT)
fit2 <- eBayes(fit2)
test.DE <- decideTests(fit2)
table(test.DE@.Data)# number of significantly differently expressed genes (down-/up-regulated)
logFC <- fit2$coefficients[,1]#logFC changes
t <- fit2$t[,1]#	numeric vector or matrix of moderated t-statistics
P.value <- fit2$p.value[,1]#numeric vector of p-values corresponding to the t-statistics
adj.P.Val <- p.adjust(fit2$p.value[,1], method="fdr")
DE <- test.DE@.Data[,1]

#Expression differences between Hot.Dry-lines and Control-lines in hot-dry conditions 
fit2 <- contrasts.fit(fit, baselineDiff.inHD)
fit2 <- eBayes(fit2)
test.DE <- decideTests(fit2)
table(test.DE@.Data)
logFC <- fit2$coefficients[,1]#logFC changes
t <- fit2$t[,1]#	numeric vector or matrix of moderated t-statistics
P.value <- fit2$p.value[,1]#numeric vector of p-values corresponding to the t-statistics
adj.P.Val <- p.adjust(fit2$p.value[,1], method="fdr")
DE <- test.DE@.Data[,1]

#Plastic response of Control-lines to hot-dry treatment 
fit2 <- contrasts.fit(fit, responseCT)
fit2 <- eBayes(fit2)
test.DE <- decideTests(fit2)
table(test.DE@.Data)
logFC <- fit2$coefficients[,1]#logFC changes
t <- fit2$t[,1]#	numeric vector or matrix of moderated t-statistics
P.value <- fit2$p.value[,1]#numeric vector of p-values corresponding to the t-statistics
adj.P.Val <- p.adjust(fit2$p.value[,1], method="fdr")
DE <- test.DE@.Data[,1]

#Plastic response of Hot.Dry-lines to hot-dry treatment 
fit2 <- contrasts.fit(fit, responseHD)
fit2 <- eBayes(fit2)
test.DE <- decideTests(fit2)
table(test.DE@.Data)
logFC <- fit2$coefficients[,1]#logFC changes
t <- fit2$t[,1]#	numeric vector or matrix of moderated t-statistics
P.value <- fit2$p.value[,1]#numeric vector of p-values corresponding to the t-statistics
adj.P.Val <- p.adjust(fit2$p.value[,1], method="fdr")
DE <- test.DE@.Data[,1]

# Differences in plasticity between Hot.Dry- and Control-lines
fit2 <- contrasts.fit(fit, Diff.responses)
fit2 <- eBayes(fit2)
test.DE <- decideTests(fit2)
table(test.DE@.Data)
logFC <- fit2$coefficients[,1]#logFC changes
t <- fit2$t[,1]#	numeric vector or matrix of moderated t-statistics
P.value <- fit2$p.value[,1]#numeric vector of p-values corresponding to the t-statistics
adj.P.Val <- p.adjust(fit2$p.value[,1], method="fdr")
DE <- test.DE@.Data[,1]

# Total difference between Control-lines in control condition
#and Hot.Dry-lines in hot-dry conditions
fit2 <- contrasts.fit(fit, total.change)
fit2 <- eBayes(fit2)
test.DE <- decideTests(fit2)
table(test.DE@.Data)
logFC <- fit2$coefficients[,1]#logFC changes
t <- fit2$t[,1]#	numeric vector or matrix of moderated t-statistics
P.value <- fit2$p.value[,1]#numeric vector of p-values corresponding to the t-statistics
adj.P.Val <- p.adjust(fit2$p.value[,1], method="fdr")
DE <- test.DE@.Data[,1]


#save(list=ls(),file="Rdata/Koch.RData")
#https://cloud.biologie.ens.fr/public.php?service=files&t=b3i8vdAevvtCtOR
