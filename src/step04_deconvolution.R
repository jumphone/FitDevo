# /home/toolkit/tools/R4.0.3/bin/R


source('https://gitee.com/jumphone/public/raw/master/fitdevo.R')
GW=readRDS(url('https://gitee.com/jumphone/public/raw/master/GW.rds'))
BGW=readRDS(url('https://gitee.com/jumphone/public/raw/master/BGW.rds'))



##########################
setwd('/home/zhangfeng/project/DEV/data/ForeskinEpidermis_GSE147482')
##########################

library(Seurat)

D1=Read10X('Donor1')
D2=Read10X('Donor2')
D3=Read10X('Donor3')
D4=Read10X('Donor4')
D5=Read10X('Donor5')
colnames(D1)=paste0('D1_',colnames(D1))
colnames(D2)=paste0('D2_',colnames(D2))
colnames(D3)=paste0('D3_',colnames(D3))
colnames(D4)=paste0('D4_',colnames(D4))
colnames(D5)=paste0('D5_',colnames(D5))


COM=.simple_combine(D1,D2)$combine
COM=.simple_combine(COM,D3)$combine
COM=.simple_combine(COM,D4)$combine
COM=.simple_combine(COM,D5)$combine

BATCH=c(rep('D1',ncol(D1)),rep('D2',ncol(D2)),rep('D3',ncol(D3)),rep('D4',ncol(D4)),rep('D5',ncol(D5)) )

pbmc=CreateSeuratObject(counts = COM, min.cells = 0, min.features = 0, project = "ALL")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc$batch=BATCH

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)


source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')

DATA=pbmc@assays$RNA@counts
BATCH=pbmc$batch

mybeer=BEER(DATA, BATCH, GNUM=30, PCNUM=50, ROUND=1, GN=2000, SEED=1, COMBAT=TRUE, RMG=NULL)   


pbmc <- mybeer$seurat
PCUSE <- mybeer$select

pbmc <- RunUMAP(object = pbmc, reduction='pca',dims = PCUSE[1:10], check_duplicates=FALSE)

COUNT=pbmc@assays$RNA@counts
RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>0),]



##################################################################
# See step02 & step03 to find ".runMethod"

SCORE_SSGW=.runMethod_SSGW(COUNT) 
SCORE_CCAT=.runMethod_CCAT(.normData(COUNT))
SCORE_CytoTRACE=.runMethod_CytoTRACE(COUNT)

pbmc$score_ssgw=SCORE_SSGW
pbmc$score_ccat=SCORE_CCAT
pbmc$score_cytotrace=SCORE_CytoTRACE

#########################################


pbmc <- FindNeighbors(pbmc, dims = PCUSE[1:10])
pbmc <- FindClusters(pbmc, resolution = 1)
DimPlot(pbmc,label=TRUE)+NoLegend()

CLST=as.numeric(as.character(pbmc$RNA_snn_res.1))

###############
# Use "FeaturePlot" to help define cell type of each cluster
FeaturePlot(pbmc, features=c('KRT5'),order=TRUE)
##################

TYPE=rep(NA,ncol(pbmc))
TYPE[which(CLST %in% c(14,16,18))]='MEL'

TYPE[which(CLST %in% c(9))]='BAS.I'
TYPE[which(CLST %in% c(7))]='BAS.II'
TYPE[which(CLST %in% c(2,19))]='BAS.III'
TYPE[which(CLST %in% c(15))]='BAS.IV'

TYPE[which(CLST %in% c(17,1,8))]='GRN'
TYPE[which(CLST %in% c(6,5,10,0,3,13,4,11,12))]='SPN'

pbmc$type=TYPE
Idents(pbmc)=factor(as.character(pbmc$type),levels=c('BAS.I','BAS.II','BAS.III','BAS.IV','SPN','GRN','MEL'))
##################################

pbmc.withoutNA=subset(pbmc, cells=colnames(pbmc)[which(!is.na(pbmc$type))])


USED=which(pbmc.withoutNA$type != 'MEL')
UDATA=pbmc.withoutNA[['RNA']]@data[,USED]

US1=pbmc.withoutNA$score_ssgw[USED]
US2=pbmc.withoutNA$score_ccat[USED]
US3=pbmc.withoutNA$score_cytotrace[USED]

BIN_ID=rep(0:9,each=ncol(UDATA)/10+1)[1:ncol(UDATA)]

BIN_ID_S1=rep(NA,ncol(UDATA))
BIN_ID_S2=rep(NA,ncol(UDATA))
BIN_ID_S3=rep(NA,ncol(UDATA))

BIN_ID_S1[order(US1)]=BIN_ID
BIN_ID_S2[order(US2)]=BIN_ID
BIN_ID_S3[order(US3)]=BIN_ID


######################################
# CIBERSORTx
###################################

USED=which(pbmc.withoutNA$type != 'MEL')
set.seed(123)
USED=sample(USED,1000)

UDATA=pbmc.withoutNA@assays$RNA@data[,USED]
RSUM=rowSums(UDATA)
UDATA=as.matrix(UDATA[which(RSUM>0),])
UDATA=apply(UDATA,2,round,3)


US1=pbmc.withoutNA$score_ssgw[USED]
US2=pbmc.withoutNA$score_ccat[USED]
US3=pbmc.withoutNA$score_cytotrace[USED]

BIN_ID=rep(0:9,each=ncol(UDATA)/10+1)[1:ncol(UDATA)]

BIN_ID_S1=rep(NA,ncol(UDATA))
BIN_ID_S2=rep(NA,ncol(UDATA))
BIN_ID_S3=rep(NA,ncol(UDATA))


BIN_ID_S1[order(US1)]=BIN_ID
BIN_ID_S2[order(US2)]=BIN_ID
BIN_ID_S3[order(US3)]=BIN_ID


OD1=UDATA
colnames(OD1)=paste0('bin',BIN_ID_S1)

OD2=UDATA
colnames(OD2)=paste0('bin',BIN_ID_S2)

OD3=UDATA
colnames(OD3)=paste0('bin',BIN_ID_S3)


write.table(OD1,file='scRNA_SSGW.txt',sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)
write.table(OD2,file='scRNA_CCAT.txt',sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)
write.table(OD3,file='scRNA_CytoTRACE.txt',sep='\t',row.names=TRUE,col.names=TRUE,quote=FALSE)




################################################



##############################################
OUT_S1=read.table('/home/zhangfeng/project/DEV/data/bulk_layer_Epidermal_GSE52651/CIBERSORTx_SSGW_array_rma.txt',sep='\t',row.names=1,header=TRUE)
OUT_S2=read.table('/home/zhangfeng/project/DEV/data/bulk_layer_Epidermal_GSE52651/CIBERSORTx_CCAT_array_rma.txt',sep='\t',row.names=1,header=TRUE)
OUT_S3=read.table('/home/zhangfeng/project/DEV/data/bulk_layer_Epidermal_GSE52651/CIBERSORTx_CytoTRACE_array_rma.txt',sep='\t',row.names=1,header=TRUE)



OUT_S1=OUT_S1[,c(1:10)]
OUT_S2=OUT_S2[,c(1:10)]
OUT_S3=OUT_S3[,c(1:10)]

OUT_S1=OUT_S1[,order(colnames(OUT_S1))]
OUT_S2=OUT_S2[,order(colnames(OUT_S2))]
OUT_S3=OUT_S3[,order(colnames(OUT_S3))]

library('ComplexHeatmap')
library('circlize')

color_fun =colorRamp2(c(-1.5,0,1.5,2 ), c('royalblue3','white','red1','gold1'))

color_fun =colorRamp2(c(-1.5,0,1.5,2 ), c('white','grey90','grey40','grey20'))

################################################
pdf('/home/zhangfeng/project/DEV/data/figure/f08_07_CIBERSORTx.pdf',width=4,height=4.5)
library('ComplexHeatmap')
library('circlize')

mat=OUT_S1
mat=as.matrix(mat)

smat=apply(mat,2,scale)
rownames(smat)=rownames(mat)
colnames(smat)=colnames(mat)

Heatmap(smat, row_title='', name="Z",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun, border = TRUE ,
          #cell_fun = function(j, i, x, y, width, height, fill) {
          #      if(smat[i, j] > .1)
          #      grid.text(sprintf("+"), x, y, gp = gpar(fontsize = 10))
               #grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 15))
              # }
	      )


################################################

mat=OUT_S2
mat=as.matrix(mat)


smat=apply(mat,2,scale)
rownames(smat)=rownames(mat)
colnames(smat)=colnames(mat)

Heatmap(smat, row_title='', name="Z",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun, border = TRUE ,
          #cell_fun = function(j, i, x, y, width, height, fill) {
           #     if(smat[i, j] > .1)
           #     grid.text(sprintf("+"), x, y, gp = gpar(fontsize = 10))
               #grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 15))
              # }
	      )

################################################

mat=OUT_S3
mat=as.matrix(mat)


smat=apply(mat,2,scale)
rownames(smat)=rownames(mat)
colnames(smat)=colnames(mat)

Heatmap(smat, row_title='', name="Z",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun, border = TRUE ,
         # cell_fun = function(j, i, x, y, width, height, fill) {
            #    if(smat[i, j] > .1)
             #   grid.text(sprintf("+"), x, y, gp = gpar(fontsize = 10))
               #grid.text(sprintf("%.1f", small_mat[i, j]), x, y, gp = gpar(fontsize = 15))
              # }
	      )


dev.off()




TIME=c(7,7,6,6,5,5,4,4,3,3,2,2,1,1,0,0)

pdf('/home/zhangfeng/project/DEV/data/figure/f08_08_CIBERSORTx_corBar.pdf',width=4,height=4)

mat=OUT_S1
mat=as.matrix(mat)
smat=t(apply(mat,1,scale))
rownames(smat)=rownames(mat)
colnames(smat)=colnames(mat)

cor(c(0:9),apply(mat,2,cor,TIME,method='pearson'))
0.711636

smat=mat

x=c(0:9)
y=apply(smat,2,cor,TIME,method='pearson')
y[which(is.na(y))]=0

plot(ylim=c(-1,1),x,y,type='p',pch='+',cex=2)
#points(x,predict(loess(y~x)),lwd=5,lty=1,type='l',col='grey80')
points(x,predict(lm(y~x)),lwd=7,lty=1,type='l',col='grey80')
points(ylim=c(-1,1),x,y,type='h',lwd=2,lty=2)
points(x,y,type='p',pch='+',cex=2)
abline(h=0,lty=2,lwd=1)

barplot(ylim=c(-1,1),apply(smat,2,cor,TIME,method='pearson'))


#####################
mat=OUT_S2
mat=as.matrix(mat)
smat=t(apply(mat,1,scale))
rownames(smat)=rownames(mat)
colnames(smat)=colnames(mat)
cor(c(0:9),apply(mat,2,cor,TIME,method='pearson'))
0.6725882

smat=mat

x=c(0:9)
y=apply(smat,2,cor,TIME,method='pearson')
y[which(is.na(y))]=0


plot(ylim=c(-1,1),x,y,type='p',pch='+',cex=2)
#points(x,predict(loess(y~x)),lwd=5,lty=1,type='l',col='grey80')
points(x,predict(lm(y~x)),lwd=7,lty=1,type='l',col='grey80')
points(ylim=c(-1,1),x,y,type='h',lwd=2,lty=2)
points(x,y,type='p',pch='+',cex=2)
abline(h=0,lty=2,lwd=1)


barplot(ylim=c(-1,1),apply(smat,2,cor,TIME,method='pearson'))



##############################
mat=OUT_S3
mat=as.matrix(mat)
smat=t(apply(mat,1,scale))
rownames(smat)=rownames(mat)
colnames(smat)=colnames(mat)
cor(c(0:9),apply(mat,2,cor,TIME,method='pearson'))
0.4148107

smat=mat

x=c(0:9)
y=apply(smat,2,cor,TIME,method='pearson')
y[which(is.na(y))]=0


plot(ylim=c(-1,1),x,y,type='p',pch='+',cex=2)
#points(x,predict(loess(y~x)),lwd=5,lty=1,type='l',col='grey80')
points(x,predict(lm(y~x)),lwd=7,lty=1,type='l',col='grey80')
points(ylim=c(-1,1),x,y,type='h',lwd=2,lty=2)
points(x,y,type='p',pch='+',cex=2)
abline(h=0,lty=2,lwd=1)


barplot(ylim=c(-1,1),apply(smat,2,cor,TIME,method='pearson'))


dev.off()








####################################################################################################
# Correlation-based deconvolution method (using HCL to conduct normalization
######################################################################################################


# /home/toolkit/tools/R4.0.3/bin/R

source('https://gitee.com/jumphone/public/raw/master/fitdevo_source.R')


GW=readRDS(file='/home/zhangfeng/project/DEV/data/RDS/GW.rds')

BGW=GW
BGW[which(GW>0)]=1
BGW[which(GW<=0)]=0
###########################

##########################
setwd('/home/zhangfeng/project/DEV/data/ForeskinEpidermis_GSE147482')
#########################

library(Seurat)

###############################
library('sva')
library('limma')
source('https://gitee.com/jumphone/Delia/raw/master/Delia.R')
source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')





pbmc=readRDS(file='/home/zhangfeng/project/DEV/data/RDS/Epidermis_seurat.rds')
pbmc.withoutNA=subset(pbmc, cells=colnames(pbmc)[which(!is.na(pbmc$type))])
USED=which(pbmc.withoutNA$type != 'MEL')
####
#without down-sampling
####
UDATA=pbmc.withoutNA@assays$RNA@data[,USED]
RSUM=rowSums(UDATA)
UDATA=as.matrix(UDATA[which(RSUM>0),])


US1=pbmc.withoutNA$score_ssgw[USED]
US2=pbmc.withoutNA$score_ccat[USED]
US3=pbmc.withoutNA$score_cytotrace[USED]

BIN_ID=rep(0:9,each=ncol(UDATA)/10+1)[1:ncol(UDATA)]

BIN_ID_S1=rep(NA,ncol(UDATA))
BIN_ID_S2=rep(NA,ncol(UDATA))
BIN_ID_S3=rep(NA,ncol(UDATA))


BIN_ID_S1[order(US1)]=BIN_ID
BIN_ID_S2[order(US2)]=BIN_ID
BIN_ID_S3[order(US3)]=BIN_ID





##########################
BULK=.readTable('/home/zhangfeng/project/DEV/data/bulk_layer_Epidermal_GSE52651/array_rma_GSE52651.txt',
      SEP='\t', UP=FALSE)
BULK=BULK[which(rowSums(BULK)>0),]

############################
HCL=readRDS('/home/zhangfeng/project/DEV/data/RDS/HCL.rds')

.logNorm10k=function(x){
    y=apply(x,2,log1p)
    csum=colSums(y)
    y=t(t(y)/csum * 10000)
    return(y)
    }

HCL=.logNorm10k(HCL)
HCL=HCL[which(rowSums(HCL)>0),]

saveRDS(HCL,'/home/zhangfeng/project/DEV/data/RDS/HCL_log10k.rds')


###########################


source('https://gitee.com/jumphone/BEER/raw/master/BEER.R')
HCL=readRDS('/home/zhangfeng/project/DEV/data/RDS/HCL_log10k.rds')

.corNormByBKG <-function( BULK, REF, BKG){
    BULK=BULK
    REF=REF
    BKG=BKG
    ######################################
    print('Combine Data ... ')
    #########################################
    COM=.simple_combine(BULK, REF)$combine
    COM=.simple_combine(COM,BKG)$combine
    BATCH=c(rep('BULK',ncol(BULK)),rep('REF',ncol(REF)),rep('BKG',ncol(BKG)))    
    #########################################
    print('Remove Batch Effect ... ')
    #########################################
    COMBAT=.combat(COM,BATCH)
    D1=COMBAT[,which(BATCH %in% c('REF','BKG'))]
    D2=COMBAT[,which(BATCH=='BULK')]
    #########################################
    print('Calculate Correlation ... ')
    #########################################
    COR=cor(D1,D2)
    SCOR=COR
    SCOR=t(apply(SCOR,1,scale))
    SCOR=apply(SCOR,2,scale)
    rownames(SCOR)=rownames(COR)
    colnames(SCOR)=colnames(COR)
    SCOR=SCOR[1:ncol(REF),]
    #######################################
    print('Finished!')
    #######################################
    OUT=t(SCOR)
    return(OUT)
    }






pdf('/home/zhangfeng/project/DEV/data/figure/f08_10_CorHCL.pdf',width=4,height=4.5)


REF=UDATA
TAG=as.character(BIN_ID_S1)


REF <- .generate_mean(REF,TAG)
REF=REF[,order(colnames(REF))]
######################
REF=REF[which(rowSums(REF)>0),]
###############################

OUT=.corNormByBKG(BULK, REF, HCL)

USED=which(pbmc.withoutNA$type != 'MEL')



TIME=c(7,7,6,6,5,5,4,4,3,3,2,2,1,1,0,0)

mat=OUT
mat=as.matrix(mat)
cor(c(0:9),apply(mat,2,cor,TIME,method='pearson'))
0.6110524

library('ComplexHeatmap')
library('circlize')
color_fun =colorRamp2(c(-1.5,0,1.5,2 ), c('white','grey90','grey40','grey20'))

smat=mat
rownames(smat)=rownames(mat)
Heatmap(smat, row_title='', name="Z",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun, border = TRUE 
	      )



##########################################################################

x=c(0:9)
y=apply(mat,2,cor,TIME,method='pearson')
y[which(is.na(y))]=0


plot(ylim=c(-1,1),x,y,type='p',pch='+',cex=2)
points(x,predict(lm(y~x)),lwd=7,lty=1,type='l',col='grey80')
points(ylim=c(-1,1),x,y,type='h',lwd=2,lty=2)
points(x,y,type='p',pch='+',cex=2)
abline(h=0,lty=2,lwd=1)

###########################################################################
###########################################################################



REF=UDATA
TAG=as.character(BIN_ID_S2)


REF <- .generate_mean(REF,TAG)
REF=REF[,order(colnames(REF))]
######################
REF=REF[which(rowSums(REF)>0),]
###############################

OUT=.corNormByBKG(BULK, REF, HCL)


TIME=c(7,7,6,6,5,5,4,4,3,3,2,2,1,1,0,0)

mat=OUT
mat=as.matrix(mat)
cor(c(0:9),apply(mat,2,cor,TIME,method='pearson'))
0.4534973

library('ComplexHeatmap')
library('circlize')
color_fun =colorRamp2(c(-1.5,0,1.5,2 ), c('white','grey90','grey40','grey20'))

smat=mat
rownames(smat)=rownames(mat)
Heatmap(smat, row_title='', name="Z",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun, border = TRUE 
	      )



##########################################################################

x=c(0:9)
y=apply(mat,2,cor,TIME,method='pearson')
y[which(is.na(y))]=0


plot(ylim=c(-1,1),x,y,type='p',pch='+',cex=2)
points(x,predict(lm(y~x)),lwd=7,lty=1,type='l',col='grey80')
points(ylim=c(-1,1),x,y,type='h',lwd=2,lty=2)
points(x,y,type='p',pch='+',cex=2)
abline(h=0,lty=2,lwd=1)




###########################################################################
###########################################################################





REF=UDATA
TAG=as.character(BIN_ID_S3)



REF <- .generate_mean(REF,TAG)
REF=REF[,order(colnames(REF))]
######################
REF=REF[which(rowSums(REF)>0),]
###############################

OUT=.corNormByBKG(BULK, REF, HCL)


TIME=c(7,7,6,6,5,5,4,4,3,3,2,2,1,1,0,0)

mat=OUT
mat=as.matrix(mat)
cor(c(0:9),apply(mat,2,cor,TIME,method='pearson'))
0.59011

library('ComplexHeatmap')
library('circlize')
color_fun =colorRamp2(c(-1.5,0,1.5,2 ), c('white','grey90','grey40','grey20'))

smat=mat
rownames(smat)=rownames(mat)
Heatmap(smat, row_title='', name="Z",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=TRUE, show_row_names=TRUE,
	      col=color_fun, border = TRUE 
	      )



##########################################################################

x=c(0:9)
y=apply(mat,2,cor,TIME,method='pearson')
y[which(is.na(y))]=0


plot(ylim=c(-1,1),x,y,type='p',pch='+',cex=2)
points(x,predict(lm(y~x)),lwd=7,lty=1,type='l',col='grey80')
points(ylim=c(-1,1),x,y,type='h',lwd=2,lty=2)
points(x,y,type='p',pch='+',cex=2)
abline(h=0,lty=2,lwd=1)


#####################################
dev.off()







