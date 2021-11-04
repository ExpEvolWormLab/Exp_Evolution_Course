#PoolSeq R package
#Link GitHub:
#  https://github.com/ThomasTaus/poolSeq/archive/refs/tags/v0.3.5.tar.gz
#Install command in R:
#  install.packages("~/Downloads/poolSeq-0.3.5.tar.gz", repos=NULL, type="source")


library('poolSeq')

## Downloaded the main genomic data file from:
#https://datadryad.org/stash/downloads/file_stream/14282

#Awk command to filter on CHR value / kept 2L only
#awk '($1=="2L")'  Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID.sync > Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID_CHR_2L.sync

# Load the file in R

#sync=read.table("Dsim_F0-F60_Q20_polymorphic_CMH_FET_blockID_CHR_2L.sync",colClasses="character")
#dim(sync)
#[1] 1068826      86
#write.table(sync[,c(1:13,64:73)],file="Dsim_F0-F60_CHR_2L_SMALL.sync",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

#mySync <- read.sync(file="data/Barghi/Dsim_F0-F60_CHR_2L_SMALL.sync", gen=rep(c(0, 60),each=10), repl=1:10)
#str(mySync)

### Large files can be downloaded here:
#https://cloud.biologie.ens.fr/public.php?service=files&t=i5VoeekmmqxbQCd

# The RData file only 100MB:
#https://cloud.biologie.ens.fr/public.php?service=files&t=G1jMzXsbyyly1PF
load("Rdata/Barghi_2L_Small.RData")


cov_F0 <- coverage(mySync,gen=0,repl=1:10)
cov_F60 <- coverage(mySync,gen=60,repl=1:10)

A0 <- af(mySync,gen=0,repl=1:10)*cov_F0
a0 <- cov_F0 - A0

At <- af(mySync,gen=60,repl=1:10)*cov_F60
at <- cov_F60 - At

hist(cov_F60)

first_CMH <- cmh.test(t(A0), t(a0), t(At), t(at), min.cov = 5, max.cov = 250, min.cnt = 1, log = FALSE)
sum(is.na(first_CMH))

table(log10(first_CMH)<(-5))
plot(-log10(first_CMH)[-log10(first_CMH)>5],pch=16,cex=.8)


