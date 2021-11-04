rm(list=ls())

library("edgeR")


data <- read.table(file="data/Huang_gene_count.txt",h=TRUE)
gene_names <- data[,1]
data <- data[,2:ncol(data)]
head(data)

keep=log2(rowSums(cpm(data,normalized.lib.sizes=TRUE))/ncol(data))>0
summary(keep)

## Number of analyzed genes
# summary(keep)
#   Mode   FALSE    TRUE 
#logical    2062   11200

## Keep the index of the spatial and temporal replicates only for analysis
names(data)
Temporal_C_index <- 27:36
Spatial_index <- 37:46

dim(data[keep,c(Temporal_C_index,Spatial_index)]) # 9957 / 20
data_ST <- data[keep,c(Temporal_C_index,Spatial_index)]
dim(data_ST)
DGEList=DGEList(data_ST)
DGEList=calcNormFactors(DGEList, method=c("TMM"))

Pop <- rep(c("Temporal","Spatial"),each=10)
Env <- rep(c("C","N"),10)
design_plasticity=model.matrix(~Pop*Env)

####### Run the GLM using edgeR

DGEList=calcNormFactors(DGEList, method=c("TMM"))
BCH_DGEList_GLM <- estimateGLMRobustDisp(DGEList, design_plasticity,verbose=TRUE,maxit=5)
BCH_fit_DGEList_GLM=glmFit(BCH_DGEList_GLM, design=design_plasticity)


## Build contrasts 

ST_inC = c(0,1,0,0)
ST_inN = c(0,1,0,1)

S_plasticity = c(0,0,1,0)
T_plasticity = c(0,0,1,1)

ST_plasticity = c(0,0,0,1)


LRT_fits_ST_inC <- glmLRT(BCH_fit_DGEList_GLM, contrast=ST_inC)	
LRT_fits_ST_inC_DE <- decideTestsDGE(LRT_fits_ST_inC, p=0.05, adjust="BH")

LRT_fits_ST_inN <- glmLRT(BCH_fit_DGEList_GLM, contrast=ST_inN)	
LRT_fits_ST_inN_DE <- decideTestsDGE(LRT_fits_ST_inN, p=0.05, adjust="BH")

table(LRT_fits_ST_inC_DE)
table(LRT_fits_ST_inN_DE)

write.table(gene_names[keep][LRT_fits_ST_inC_DE==c(-1)],file="txt/Huang_ST_inC_dec.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(gene_names[keep],file="txt/Huang_all_analysed.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)

#http://cbl-gorilla.cs.technion.ac.il/


LRT_fits_T_plasticity <- glmLRT(BCH_fit_DGEList_GLM, contrast=T_plasticity)	
LRT_fits_T_plasticity_DE <- decideTestsDGE(LRT_fits_T_plasticity, p=0.05, adjust="BH")

table(LRT_fits_T_plasticity_DE)

write.table(gene_names[keep][LRT_fits_T_plasticity_DE==c(-1)],file="txt/Huang_T_plasticity_dec.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)
write.table(gene_names[keep][LRT_fits_T_plasticity_DE==c(1)],file="txt/Huang_T_plasticity_inc.txt",row.names=FALSE,col.names=FALSE,quote=FALSE)



