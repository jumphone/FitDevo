# /home/toolkit/tools/R4.0.3/bin/R


##########################
setwd('/home/zhangfeng/project/DEV/data/breastCancer_NatureGenetics/Wu_etal_2021_BRCA_scRNASeq')
##################



library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=TRUE)


GW=readRDS(file='/home/zhangfeng/project/DEV/data/RDS/GW.rds')

BGW=GW
BGW[which(GW>0)]=1
BGW[which(GW<=0)]=0
###########################

.getBin<-function(x,n=10){
    if(length(x) %% n !=0){
        BIN=rep(c(0:(n-1)),each=(length(x)/n+1) )
        }else{
        BIN=rep(c(0:(n-1)),each=(length(x)/n) )
        }
    y=x
    y[order(x)]=BIN[1:length(x)]
    return(y)
    }



library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratDisk)

DATA=Read10X('/home/zhangfeng/project/DEV/data/breastCancer_NatureGenetics/Wu_etal_2021_BRCA_scRNASeq',gene.column =1)
META=read.table('/home/zhangfeng/project/DEV/data/breastCancer_NatureGenetics/Wu_etal_2021_BRCA_scRNASeq/metadata.csv',sep=',',header=TRUE,row.names=1)



pbmc=CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = "ALL",meta.data=META)
#100064

TYPE1=pbmc$celltype_major
TYPE2=pbmc$celltype_minor

T_INDEX=which(TYPE1=='Cancer Epithelial')
pbmc.tumor=subset(pbmc,cells=colnames(pbmc)[T_INDEX])

COUNT=pbmc.tumor[['RNA']]@counts

RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>0),]
#################
SCORE_SSGW=.runMethod_SSGW(COUNT)
SCORE_CCAT=.runMethod_CCAT(.normData(COUNT))
SCORE_CytoTRACE=.runMethod_CytoTRACE(COUNT)


pbmc.tumor$score_ssgw=SCORE_SSGW
pbmc.tumor$score_ccat=SCORE_CCAT
pbmc.tumor$score_cytotrace=SCORE_CytoTRACE

   
bin_ssgw=.getBin(pbmc.tumor$score_ssgw,10)
bin_ccat=.getBin(pbmc.tumor$score_ccat,10)
bin_cytotrace=.getBin(pbmc.tumor$score_cytotrace,10)


#################################
run_deseq2<-function(data,tag,level){
    library(DESeq2)
    mydata=data
    type<-factor(tag,levels = level)
    database<-as.matrix(mydata)
    coldata<-data.frame(row.names = colnames(database),type)
    dds<-DESeqDataSetFromMatrix(database,coldata,design = ~type)
    ###################################################
    CPU=10
    BPPARAM = BatchJobsParam(CPU)
    dds<-DESeq(dds,parallel=TRUE,BPPARAM=BPPARAM)
    ###################################################
    #dds<-DESeq(dds)
    res = results(dds, contrast=c("type", c(level)))
    res=as.matrix(res)
    return(res)
    }
########################################################
COUNT=pbmc.tumor[['RNA']]@counts
RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>100),]
###################







#############################
# Use three different seeds (123,132,321) to generate pseudo-bulk and identify signatures
##############################

SEED=123



###########################################################################################################
BIN=bin_ssgw
###########################################################################################################
Idents(pbmc.tumor)=factor(BIN,levels=sort(unique(BIN)))

UBIN=sort(unique(BIN))

FOLD=rep(10,length(UBIN))
names(FOLD)=UBIN
mybeer.agg=BEER.AGG(COUNT, BIN, FOLD, PCNUM=50, GN=2000, SEED=SEED)

COUNT.AGG=mybeer.agg$data.agg
BIN.AGG=mybeer.agg$data.agg.batch


OUT=list()
i=1
while(i<=length(UBIN)){
    this_bin=UBIN[i]
    print(this_bin)
    this_index=which(BIN.AGG==this_bin)
    all_index=c(1:length(BIN.AGG))
    used_bkg_index=all_index[which(!all_index %in% this_index)]
    this_count=COUNT.AGG[,c(this_index,used_bkg_index)]
    colnames(this_count)[c(1:length(this_index))]=paste0('target.', colnames(this_count)[c(1:length(this_index))])
    this_count=as.matrix(this_count[which(rowSums(this_count)>0),])
    this_tag=c(rep('target',length(this_index)),rep('bkg',length(used_bkg_index)))
    this_level=c('target','bkg')
    MK=run_deseq2(this_count,this_tag,this_level)
    UMK=MK
    UMK=UMK[order(-UMK[,2]),]
    print(head(UMK,n=10))
    OUT[[i]]=UMK
    i=i+1
    }


saveRDS(OUT, file=paste0('deseq2_ssgw_agg',SEED,'.rds'))
saveRDS(COUNT.AGG, file=paste0('deseq2_ssgw_agg',SEED,'.count.rds'))
saveRDS(BIN.AGG, file=paste0('deseq2_ssgw_agg',SEED,'.bin.rds'))

####################################################################################################################################################







###########################################################################################################
BIN=bin_ccat
###########################################################################################################
Idents(pbmc.tumor)=factor(BIN,levels=sort(unique(BIN)))

UBIN=sort(unique(BIN))

FOLD=rep(10,length(UBIN))
names(FOLD)=UBIN
mybeer.agg=BEER.AGG(COUNT, BIN, FOLD, PCNUM=50, GN=2000, SEED=SEED)

COUNT.AGG=mybeer.agg$data.agg
BIN.AGG=mybeer.agg$data.agg.batch


OUT=list()
i=1
while(i<=length(UBIN)){
    this_bin=UBIN[i]
    print(this_bin)
    this_index=which(BIN.AGG==this_bin)
    all_index=c(1:length(BIN.AGG))
    used_bkg_index=all_index[which(!all_index %in% this_index)]
    this_count=COUNT.AGG[,c(this_index,used_bkg_index)]
    colnames(this_count)[c(1:length(this_index))]=paste0('target.', colnames(this_count)[c(1:length(this_index))])
    this_count=as.matrix(this_count[which(rowSums(this_count)>0),])
    this_tag=c(rep('target',length(this_index)),rep('bkg',length(used_bkg_index)))
    this_level=c('target','bkg')
    MK=run_deseq2(this_count,this_tag,this_level)
    UMK=MK
    UMK=UMK[order(-UMK[,2]),]
    print(head(UMK,n=10))
    OUT[[i]]=UMK
    i=i+1
    }


saveRDS(OUT, file=paste0('deseq2_ccat_agg',SEED,'.rds'))
saveRDS(COUNT.AGG, file=paste0('deseq2_ccat_agg',SEED,'.count.rds'))
saveRDS(BIN.AGG, file=paste0('deseq2_ccat_agg',SEED,'.bin.rds'))

####################################################################################################################################################



###########################################################################################################
BIN=bin_cytotrace
###########################################################################################################
Idents(pbmc.tumor)=factor(BIN,levels=sort(unique(BIN)))

UBIN=sort(unique(BIN))

FOLD=rep(10,length(UBIN))
names(FOLD)=UBIN
mybeer.agg=BEER.AGG(COUNT, BIN, FOLD, PCNUM=50, GN=2000, SEED=SEED)

COUNT.AGG=mybeer.agg$data.agg
BIN.AGG=mybeer.agg$data.agg.batch


OUT=list()
i=1
while(i<=length(UBIN)){
    this_bin=UBIN[i]
    print(this_bin)
    this_index=which(BIN.AGG==this_bin)
    all_index=c(1:length(BIN.AGG))
    used_bkg_index=all_index[which(!all_index %in% this_index)]
    this_count=COUNT.AGG[,c(this_index,used_bkg_index)]
    colnames(this_count)[c(1:length(this_index))]=paste0('target.', colnames(this_count)[c(1:length(this_index))])
    this_count=as.matrix(this_count[which(rowSums(this_count)>0),])
    this_tag=c(rep('target',length(this_index)),rep('bkg',length(used_bkg_index)))
    this_level=c('target','bkg')
    MK=run_deseq2(this_count,this_tag,this_level)
    UMK=MK
    UMK=UMK[order(-UMK[,2]),]
    print(head(UMK,n=10))
    OUT[[i]]=UMK
    i=i+1
    }


saveRDS(OUT, file=paste0('deseq2_cytotrace_agg',SEED,'.rds'))
saveRDS(COUNT.AGG, file=paste0('deseq2_cytotrace_agg',SEED,'.count.rds'))
saveRDS(BIN.AGG, file=paste0('deseq2_cytotrace_agg',SEED,'.bin.rds'))

####################################################################################################################################################



SEED=132



###########################################################################################################
BIN=bin_ssgw
###########################################################################################################
Idents(pbmc.tumor)=factor(BIN,levels=sort(unique(BIN)))

UBIN=sort(unique(BIN))

FOLD=rep(10,length(UBIN))
names(FOLD)=UBIN
mybeer.agg=BEER.AGG(COUNT, BIN, FOLD, PCNUM=50, GN=2000, SEED=SEED)

COUNT.AGG=mybeer.agg$data.agg
BIN.AGG=mybeer.agg$data.agg.batch


OUT=list()
i=1
while(i<=length(UBIN)){
    this_bin=UBIN[i]
    print(this_bin)
    this_index=which(BIN.AGG==this_bin)
    all_index=c(1:length(BIN.AGG))
    used_bkg_index=all_index[which(!all_index %in% this_index)]
    this_count=COUNT.AGG[,c(this_index,used_bkg_index)]
    colnames(this_count)[c(1:length(this_index))]=paste0('target.', colnames(this_count)[c(1:length(this_index))])
    this_count=as.matrix(this_count[which(rowSums(this_count)>0),])
    this_tag=c(rep('target',length(this_index)),rep('bkg',length(used_bkg_index)))
    this_level=c('target','bkg')
    MK=run_deseq2(this_count,this_tag,this_level)
    UMK=MK
    UMK=UMK[order(-UMK[,2]),]
    print(head(UMK,n=10))
    OUT[[i]]=UMK
    i=i+1
    }


saveRDS(OUT, file=paste0('deseq2_ssgw_agg',SEED,'.rds'))
saveRDS(COUNT.AGG, file=paste0('deseq2_ssgw_agg',SEED,'.count.rds'))
saveRDS(BIN.AGG, file=paste0('deseq2_ssgw_agg',SEED,'.bin.rds'))

####################################################################################################################################################







###########################################################################################################
BIN=bin_ccat
###########################################################################################################
Idents(pbmc.tumor)=factor(BIN,levels=sort(unique(BIN)))

UBIN=sort(unique(BIN))

FOLD=rep(10,length(UBIN))
names(FOLD)=UBIN
mybeer.agg=BEER.AGG(COUNT, BIN, FOLD, PCNUM=50, GN=2000, SEED=SEED)

COUNT.AGG=mybeer.agg$data.agg
BIN.AGG=mybeer.agg$data.agg.batch


OUT=list()
i=1
while(i<=length(UBIN)){
    this_bin=UBIN[i]
    print(this_bin)
    this_index=which(BIN.AGG==this_bin)
    all_index=c(1:length(BIN.AGG))
    used_bkg_index=all_index[which(!all_index %in% this_index)]
    this_count=COUNT.AGG[,c(this_index,used_bkg_index)]
    colnames(this_count)[c(1:length(this_index))]=paste0('target.', colnames(this_count)[c(1:length(this_index))])
    this_count=as.matrix(this_count[which(rowSums(this_count)>0),])
    this_tag=c(rep('target',length(this_index)),rep('bkg',length(used_bkg_index)))
    this_level=c('target','bkg')
    MK=run_deseq2(this_count,this_tag,this_level)
    UMK=MK
    UMK=UMK[order(-UMK[,2]),]
    print(head(UMK,n=10))
    OUT[[i]]=UMK
    i=i+1
    }


saveRDS(OUT, file=paste0('deseq2_ccat_agg',SEED,'.rds'))
saveRDS(COUNT.AGG, file=paste0('deseq2_ccat_agg',SEED,'.count.rds'))
saveRDS(BIN.AGG, file=paste0('deseq2_ccat_agg',SEED,'.bin.rds'))

####################################################################################################################################################



###########################################################################################################
BIN=bin_cytotrace
###########################################################################################################
Idents(pbmc.tumor)=factor(BIN,levels=sort(unique(BIN)))

UBIN=sort(unique(BIN))

FOLD=rep(10,length(UBIN))
names(FOLD)=UBIN
mybeer.agg=BEER.AGG(COUNT, BIN, FOLD, PCNUM=50, GN=2000, SEED=SEED)

COUNT.AGG=mybeer.agg$data.agg
BIN.AGG=mybeer.agg$data.agg.batch


OUT=list()
i=1
while(i<=length(UBIN)){
    this_bin=UBIN[i]
    print(this_bin)
    this_index=which(BIN.AGG==this_bin)
    all_index=c(1:length(BIN.AGG))
    used_bkg_index=all_index[which(!all_index %in% this_index)]
    this_count=COUNT.AGG[,c(this_index,used_bkg_index)]
    colnames(this_count)[c(1:length(this_index))]=paste0('target.', colnames(this_count)[c(1:length(this_index))])
    this_count=as.matrix(this_count[which(rowSums(this_count)>0),])
    this_tag=c(rep('target',length(this_index)),rep('bkg',length(used_bkg_index)))
    this_level=c('target','bkg')
    MK=run_deseq2(this_count,this_tag,this_level)
    UMK=MK
    UMK=UMK[order(-UMK[,2]),]
    print(head(UMK,n=10))
    OUT[[i]]=UMK
    i=i+1
    }


saveRDS(OUT, file=paste0('deseq2_cytotrace_agg',SEED,'.rds'))
saveRDS(COUNT.AGG, file=paste0('deseq2_cytotrace_agg',SEED,'.count.rds'))
saveRDS(BIN.AGG, file=paste0('deseq2_cytotrace_agg',SEED,'.bin.rds'))

####################################################################################################################################################



SEED=321



###########################################################################################################
BIN=bin_ssgw
###########################################################################################################
Idents(pbmc.tumor)=factor(BIN,levels=sort(unique(BIN)))

UBIN=sort(unique(BIN))

FOLD=rep(10,length(UBIN))
names(FOLD)=UBIN
mybeer.agg=BEER.AGG(COUNT, BIN, FOLD, PCNUM=50, GN=2000, SEED=SEED)

COUNT.AGG=mybeer.agg$data.agg
BIN.AGG=mybeer.agg$data.agg.batch


OUT=list()
i=1
while(i<=length(UBIN)){
    this_bin=UBIN[i]
    print(this_bin)
    this_index=which(BIN.AGG==this_bin)
    all_index=c(1:length(BIN.AGG))
    used_bkg_index=all_index[which(!all_index %in% this_index)]
    this_count=COUNT.AGG[,c(this_index,used_bkg_index)]
    colnames(this_count)[c(1:length(this_index))]=paste0('target.', colnames(this_count)[c(1:length(this_index))])
    this_count=as.matrix(this_count[which(rowSums(this_count)>0),])
    this_tag=c(rep('target',length(this_index)),rep('bkg',length(used_bkg_index)))
    this_level=c('target','bkg')
    MK=run_deseq2(this_count,this_tag,this_level)
    UMK=MK
    UMK=UMK[order(-UMK[,2]),]
    print(head(UMK,n=10))
    OUT[[i]]=UMK
    i=i+1
    }


saveRDS(OUT, file=paste0('deseq2_ssgw_agg',SEED,'.rds'))
saveRDS(COUNT.AGG, file=paste0('deseq2_ssgw_agg',SEED,'.count.rds'))
saveRDS(BIN.AGG, file=paste0('deseq2_ssgw_agg',SEED,'.bin.rds'))

####################################################################################################################################################







###########################################################################################################
BIN=bin_ccat
###########################################################################################################
Idents(pbmc.tumor)=factor(BIN,levels=sort(unique(BIN)))

UBIN=sort(unique(BIN))

FOLD=rep(10,length(UBIN))
names(FOLD)=UBIN
mybeer.agg=BEER.AGG(COUNT, BIN, FOLD, PCNUM=50, GN=2000, SEED=SEED)

COUNT.AGG=mybeer.agg$data.agg
BIN.AGG=mybeer.agg$data.agg.batch


OUT=list()
i=1
while(i<=length(UBIN)){
    this_bin=UBIN[i]
    print(this_bin)
    this_index=which(BIN.AGG==this_bin)
    all_index=c(1:length(BIN.AGG))
    used_bkg_index=all_index[which(!all_index %in% this_index)]
    this_count=COUNT.AGG[,c(this_index,used_bkg_index)]
    colnames(this_count)[c(1:length(this_index))]=paste0('target.', colnames(this_count)[c(1:length(this_index))])
    this_count=as.matrix(this_count[which(rowSums(this_count)>0),])
    this_tag=c(rep('target',length(this_index)),rep('bkg',length(used_bkg_index)))
    this_level=c('target','bkg')
    MK=run_deseq2(this_count,this_tag,this_level)
    UMK=MK
    UMK=UMK[order(-UMK[,2]),]
    print(head(UMK,n=10))
    OUT[[i]]=UMK
    i=i+1
    }


saveRDS(OUT, file=paste0('deseq2_ccat_agg',SEED,'.rds'))
saveRDS(COUNT.AGG, file=paste0('deseq2_ccat_agg',SEED,'.count.rds'))
saveRDS(BIN.AGG, file=paste0('deseq2_ccat_agg',SEED,'.bin.rds'))

####################################################################################################################################################



###########################################################################################################
BIN=bin_cytotrace
###########################################################################################################
Idents(pbmc.tumor)=factor(BIN,levels=sort(unique(BIN)))

UBIN=sort(unique(BIN))

FOLD=rep(10,length(UBIN))
names(FOLD)=UBIN
mybeer.agg=BEER.AGG(COUNT, BIN, FOLD, PCNUM=50, GN=2000, SEED=SEED)

COUNT.AGG=mybeer.agg$data.agg
BIN.AGG=mybeer.agg$data.agg.batch


OUT=list()
i=1
while(i<=length(UBIN)){
    this_bin=UBIN[i]
    print(this_bin)
    this_index=which(BIN.AGG==this_bin)
    all_index=c(1:length(BIN.AGG))
    used_bkg_index=all_index[which(!all_index %in% this_index)]
    this_count=COUNT.AGG[,c(this_index,used_bkg_index)]
    colnames(this_count)[c(1:length(this_index))]=paste0('target.', colnames(this_count)[c(1:length(this_index))])
    this_count=as.matrix(this_count[which(rowSums(this_count)>0),])
    this_tag=c(rep('target',length(this_index)),rep('bkg',length(used_bkg_index)))
    this_level=c('target','bkg')
    MK=run_deseq2(this_count,this_tag,this_level)
    UMK=MK
    UMK=UMK[order(-UMK[,2]),]
    print(head(UMK,n=10))
    OUT[[i]]=UMK
    i=i+1
    }


saveRDS(OUT, file=paste0('deseq2_cytotrace_agg',SEED,'.rds'))
saveRDS(COUNT.AGG, file=paste0('deseq2_cytotrace_agg',SEED,'.count.rds'))
saveRDS(BIN.AGG, file=paste0('deseq2_cytotrace_agg',SEED,'.bin.rds'))

####################################################################################################################################################







#########################


# /home/toolkit/tools/R4.0.3/bin/R

source('https://gitee.com/jumphone/public/raw/master/fitdevo_source.R')

library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=TRUE)


GW=readRDS(file='/home/zhangfeng/project/DEV/data/RDS/GW.rds')

BGW=GW
BGW[which(GW>0)]=1
BGW[which(GW<=0)]=0
###########################

.getBin<-function(x,n=10){
    if(length(x) %% n !=0){
        BIN=rep(c(0:(n-1)),each=(length(x)/n+1) )
        }else{
        BIN=rep(c(0:(n-1)),each=(length(x)/n) )
        }
    y=x
    y[order(x)]=BIN[1:length(x)]
    return(y)
    }



##########################
setwd('/home/zhangfeng/project/DEV/data/breastCancer_NatureGenetics/Wu_etal_2021_BRCA_scRNASeq')
##########################

library(dplyr)
library(Seurat)
library(patchwork)


pbmc.tumor=readRDS('/home/zhangfeng/project/DEV/data/breastCancer_NatureGenetics/Wu_etal_2021_BRCA_scRNASeq/pbmc.tumor.rds')



library(stringr)
##########################
BULK=.readTable('/home/zhangfeng/project/DEV/data/breastCancer_NatureGenetics/Wu_etal_2021_BRCA_scRNASeq/TCGA/TCGA-BRCA.htseq_fpkm.tsv.gene.tsv',
      SEP='\t', UP=FALSE)
BULK=BULK[which(rowSums(BULK)>0),]

SUR = read.table('/home/zhangfeng/project/DEV/data/breastCancer_NatureGenetics/Wu_etal_2021_BRCA_scRNASeq/TCGA/TCGA-BRCA.survival.tsv',sep='\t',
      row.names=1,header=TRUE) 

rownames(SUR)=str_replace_all(rownames(SUR),'-','\\.')

PHE = read.csv('/home/zhangfeng/project/DEV/data/breastCancer_NatureGenetics/Wu_etal_2021_BRCA_scRNASeq/TCGA/TCGA-BRCA.GDC_phenotype.tsv',sep='\t',
      row.names=1,header=TRUE) 
rownames(PHE)=str_replace_all(rownames(PHE),'-','\\.')

#################


CCC=.simple_combine(SUR,PHE)
DD1=as.data.frame(CCC$exp_sc_mat1)
DD2=as.data.frame(CCC$exp_sc_mat2)

TTT=rep(0,length(DD2$sample_type.samples))
TTT[which(DD2$sample_type.samples=='Primary Tumor')]=1


SSS=cbind(DD1$OS ,DD1$OS.time, DD2$age_at_diagnosis.diagnoses,TTT )
rownames(SSS)=rownames(DD1)
colnames(SSS)=c('OS','OS.time','AGE','Type')
SSS=SSS[which(!is.na(SSS[,3])),]
SSS=SSS[which(SSS[,4]==1),]


COM=.simple_combine(t(BULK),SSS)
D1=COM$exp_sc_mat1
D2=COM$exp_sc_mat2


AGE=D2[,3]


LW=quantile(AGE,0.25)
UP=quantile(AGE,0.75)

USED=which(AGE>LW & AGE< UP)



DEN=density(AGE)
COL=rep('grey90',length(DEN$x))
COL[which(DEN$x>LW & DEN$x<UP)]='grey80'
plot(DEN$x,DEN$y,col=COL,type='h')
points(DEN$x,DEN$y,type='l',lwd=2)
abline(v=LW,lwd=2,lty=2)
abline(v=UP,lwd=2,lty=2)

LW/365
UP/365


UD1=D1[USED,]
UD2=D2[USED,]





HR_LIST=readRDS('/home/zhangfeng/project/DEV/data/RDS/TCGA_BRCA_HR_LIST.rds')



###################
pdf('/home/zhangfeng/project/DEV/data/figure/f10_01_TCGA_HR.pdf',width=4,height=4)




##################################
OUT1=readRDS(file='deseq2_ssgw_agg123.rds')
OUT2=readRDS(file='deseq2_ssgw_agg321.rds')
OUT3=readRDS(file='deseq2_ssgw_agg132.rds')
###################################


FCUT=log(2,2)
PCUT=0.05
MCUT=0.5

SIG_LIST=list()
N1=c()
N2=c()
N3=c()
NI=c()

i=1
while(i<=length(OUT1)){
    this_out1=OUT1[[i]]
    this_out2=OUT2[[i]]
    this_out3=OUT3[[i]]
    ##############################################
    this_out1=this_out1[which(this_out1[,1]>MCUT),]
    this_out2=this_out2[which(this_out2[,1]>MCUT),]
    this_out3=this_out3[which(this_out3[,1]>MCUT),]
    #################################################
    this_out1=this_out1[which(this_out1[,2]>FCUT),]
    this_out2=this_out2[which(this_out2[,2]>FCUT),]
    this_out3=this_out3[which(this_out3[,2]>FCUT),]
    #################################################
    this_out1=this_out1[which(this_out1[,6]<PCUT),]
    this_out2=this_out2[which(this_out2[,6]<PCUT),]
    this_out3=this_out3[which(this_out3[,6]<PCUT),]
    ##################################################
    this_out1=this_out1[order(this_out1[,5]),]
    this_out2=this_out2[order(this_out2[,5]),]
    this_out3=this_out3[order(this_out3[,5]),]
    ##################################################
    this_gene1=rownames(this_out1)
    this_gene2=rownames(this_out2)
    this_gene3=rownames(this_out3)
    ##################################################
    this_gene=intersect(intersect(this_gene1,this_gene2),this_gene3)
    #####################################################
    N1=c(N1,length(this_gene1))
    N2=c(N2,length(this_gene2))
    N3=c(N3,length(this_gene3))
    NI=c(NI,length(this_gene))
    ###################################################
    SIG_LIST[[i]]=this_gene
    i=i+1
    }
###########################
    

#####################
XXX=c()
i=1
while(i<=length(OUT)){
    this_gene=SIG_LIST[[i]]
    this_xxx=apply(UD1[,which(colnames(UD1) %in% this_gene)],1,mean)
    XXX=cbind(XXX,this_xxx)
    i=i+1
    }


library(survival)
library(survminer)
library(survcomp)

INPUT=list()
INPUT$status=UD2[,1]
INPUT$time=UD2[,2]
INPUT$age=UD2[,3]



HH=c()
i=1
while(i<=10){
    SX=XXX[,i]
    INPUT.HR=INPUT
    INPUT.HR$cluster=rep(0,length(INPUT$time))
    INPUT.HR$cluster[which(SX>median(SX))]=1
    HR=hazard.ratio(INPUT.HR$cluster, INPUT.HR$time,INPUT.HR$status)
    HH=c(HH,HR$hazard.ratio)
    i=i+1
    }
  
cor(HH,1:10)
0.7165967
HH[10]
1.616086

plot(c(0:9),HH,type='p',pch=21,ylim=c(0.2,1.8),lwd=0,cex=0)
abline(lm(HH~c(0:9)),lwd=7,,col='grey80')
points(x=c(0:9),y=HH,type='p',pch='+',cex=2)
abline(h=1,lty=2,lwd=1)



#################################################################




##################################
OUT1=readRDS(file='deseq2_ccat_agg123.rds')
OUT2=readRDS(file='deseq2_ccat_agg321.rds')
OUT3=readRDS(file='deseq2_ccat_agg132.rds')
###################################



FCUT=log(2,2)
PCUT=0.05
MCUT=0.5

SIG_LIST=list()
N1=c()
N2=c()
N3=c()
NI=c()

i=1
while(i<=length(OUT1)){
    this_out1=OUT1[[i]]
    this_out2=OUT2[[i]]
    this_out3=OUT3[[i]]
    ##############################################
    this_out1=this_out1[which(this_out1[,1]>MCUT),]
    this_out2=this_out2[which(this_out2[,1]>MCUT),]
    this_out3=this_out3[which(this_out3[,1]>MCUT),]
    #################################################
    this_out1=this_out1[which(this_out1[,2]>FCUT),]
    this_out2=this_out2[which(this_out2[,2]>FCUT),]
    this_out3=this_out3[which(this_out3[,2]>FCUT),]
    #################################################
    this_out1=this_out1[which(this_out1[,6]<PCUT),]
    this_out2=this_out2[which(this_out2[,6]<PCUT),]
    this_out3=this_out3[which(this_out3[,6]<PCUT),]
    ##################################################
    this_out1=this_out1[order(this_out1[,5]),]
    this_out2=this_out2[order(this_out2[,5]),]
    this_out3=this_out3[order(this_out3[,5]),]
    ##################################################
    this_gene1=rownames(this_out1)
    this_gene2=rownames(this_out2)
    this_gene3=rownames(this_out3)
    ##################################################
    this_gene=intersect(intersect(this_gene1,this_gene2),this_gene3)
    #####################################################
    N1=c(N1,length(this_gene1))
    N2=c(N2,length(this_gene2))
    N3=c(N3,length(this_gene3))
    NI=c(NI,length(this_gene))
    ###################################################
    SIG_LIST[[i]]=this_gene
    i=i+1
    }
###########################
    

#####################
XXX=c()
i=1
while(i<=length(OUT1)){
    this_gene=SIG_LIST[[i]]
    this_xxx=apply(UD1[,which(colnames(UD1) %in% this_gene)],1,mean)
    XXX=cbind(XXX,this_xxx)
    i=i+1
    }


library(survival)
library(survminer)
library(survcomp)

INPUT=list()
INPUT$status=UD2[,1]
INPUT$time=UD2[,2]
INPUT$age=UD2[,3]



HH=c()
i=1
while(i<=10){
    SX=XXX[,i]
    INPUT.HR=INPUT
    INPUT.HR$cluster=rep(0,length(INPUT$time))
    INPUT.HR$cluster[which(SX>median(SX))]=1
    HR=hazard.ratio(INPUT.HR$cluster, INPUT.HR$time,INPUT.HR$status)
    HH=c(HH,HR$hazard.ratio)
    i=i+1
    }
  
cor(HH,1:10)
#0.6926282
HH[10]
#1.389603

plot(c(0:9),HH,type='p',pch=21,ylim=c(0.2,1.8),lwd=0,cex=0)
abline(lm(HH~c(0:9)),lwd=7,,col='grey80')
points(x=c(0:9),y=HH,type='p',pch='+',cex=2)
abline(h=1,lty=2,lwd=1)





##################################
OUT1=readRDS(file='deseq2_cytotrace_agg123.rds')
OUT2=readRDS(file='deseq2_cytotrace_agg321.rds')
OUT3=readRDS(file='deseq2_cytotrace_agg132.rds')
###################################



FCUT=log(2,2)
PCUT=0.05
MCUT=0.5

SIG_LIST=list()
N1=c()
N2=c()
N3=c()
NI=c()

i=1
while(i<=length(OUT1)){
    this_out1=OUT1[[i]]
    this_out2=OUT2[[i]]
    this_out3=OUT3[[i]]
    ##############################################
    this_out1=this_out1[which(this_out1[,1]>MCUT),]
    this_out2=this_out2[which(this_out2[,1]>MCUT),]
    this_out3=this_out3[which(this_out3[,1]>MCUT),]
    #################################################
    this_out1=this_out1[which(this_out1[,2]>FCUT),]
    this_out2=this_out2[which(this_out2[,2]>FCUT),]
    this_out3=this_out3[which(this_out3[,2]>FCUT),]
    #################################################
    this_out1=this_out1[which(this_out1[,6]<PCUT),]
    this_out2=this_out2[which(this_out2[,6]<PCUT),]
    this_out3=this_out3[which(this_out3[,6]<PCUT),]
    ##################################################
    this_out1=this_out1[order(this_out1[,5]),]
    this_out2=this_out2[order(this_out2[,5]),]
    this_out3=this_out3[order(this_out3[,5]),]
    ##################################################
    this_gene1=rownames(this_out1)
    this_gene2=rownames(this_out2)
    this_gene3=rownames(this_out3)
    ##################################################
    this_gene=intersect(intersect(this_gene1,this_gene2),this_gene3)
    #####################################################
    N1=c(N1,length(this_gene1))
    N2=c(N2,length(this_gene2))
    N3=c(N3,length(this_gene3))
    NI=c(NI,length(this_gene))
    ###################################################
    SIG_LIST[[i]]=this_gene
    i=i+1
    }
###########################
    

#####################
XXX=c()
i=1
while(i<=length(OUT1)){
    this_gene=SIG_LIST[[i]]
    this_xxx=apply(UD1[,which(colnames(UD1) %in% this_gene)],1,mean)
    XXX=cbind(XXX,this_xxx)
    i=i+1
    }


library(survival)
library(survminer)
library(survcomp)

INPUT=list()
INPUT$status=UD2[,1]
INPUT$time=UD2[,2]
INPUT$age=UD2[,3]



HH=c()
i=1
while(i<=10){
    SX=XXX[,i]
    INPUT.HR=INPUT
    INPUT.HR$cluster=rep(0,length(INPUT$time))
    INPUT.HR$cluster[which(SX>median(SX))]=1
    HR=hazard.ratio(INPUT.HR$cluster, INPUT.HR$time,INPUT.HR$status)
    HH=c(HH,HR$hazard.ratio)
    i=i+1
    }
  
cor(HH,1:10)
#0.6376097
HH[10]
#1.22753

plot(c(0:9),HH,type='p',pch=21,ylim=c(0.2,1.8),lwd=0,cex=0)
abline(lm(HH~c(0:9)),lwd=7,,col='grey80')
points(x=c(0:9),y=HH,type='p',pch='+',cex=2)
abline(h=1,lty=2,lwd=1)



dev.off()




















