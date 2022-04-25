

# /home/toolkit/tools/R4.0.3/bin/R

library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=TRUE)

GW=readRDS(file='/home/zhangfeng/project/DEV/data/RDS/GW.rds')

BGW=GW
BGW[which(GW>0)]=1
BGW[which(GW<=0)]=0

##########################
setwd('/home/zhangfeng/project/DEV/data/chickenHeart_GSE149457')
##########################



D4=Read10X_h5('GSM4502482_chicken_heart_spatial_RNAseq_D4_filtered_feature_bc_matrix.h5')
D7=Read10X_h5('GSM4502483_chicken_heart_spatial_RNAseq_D7_filtered_feature_bc_matrix.h5')
D10=Read10X_h5('GSM4502484_chicken_heart_spatial_RNAseq_D10_filtered_feature_bc_matrix.h5')
D14=Read10X_h5('GSM4502485_chicken_heart_spatial_RNAseq_D14_filtered_feature_bc_matrix.h5')

D4=as.matrix(D4)
D7=as.matrix(D7)
D10=as.matrix(D10)
D14=as.matrix(D14)

COM=.simple_combine(D4,D7)$combine
COM=.simple_combine(COM,D10)$combine
COM=.simple_combine(COM,D14)$combine
TIME=c(rep(4,ncol(D4)),rep(3,ncol(D7)),rep(2,ncol(D10)),rep(1,ncol(D14)) )




COUNT=COM
DATA=.normData(COUNT)


SCORE_SSGW=.runMethod_SSGW(DATA, NORM=FALSE)
SCORE_CCAT=.runMethod_CCAT(DATA)
SCORE_CytoTRACE =.runMethod_CytoTRACE(COUNT)


boxplot(SCORE_SSGW~TIME)
boxplot(SCORE_CCAT~TIME)
boxplot(SCORE_CytoTRACE~TIME)




cor(SCORE_SSGW,TIME,method='spearman')
# 0.4318146
cor(SCORE_CCAT,TIME,method='spearman')
# 0.3028755
cor(SCORE_CytoTRACE,TIME,method='spearman')
# 0.122123




###################################################

this.mat=D14
this.image=Read10X_Image('/home/zhangfeng/project/DEV/data/chickenHeart_GSE149457/chicken_heart-master/data/chicken_heart_spatial_RNAseq_processed/D14')
DefaultAssay(object = this.image) <- 'Spatial'

this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial')
this.obj[['image']] <- this.image

this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)

plot1 <- VlnPlot(this.obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(this.obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)


this.obj$score_ssgw=SCORE_SSGW[which(TIME==1)]
this.obj$score_ccat=SCORE_CCAT[which(TIME==1)]
this.obj$score_cytotrace=SCORE_CytoTRACE[which(TIME==1)]


SpatialFeaturePlot(this.obj, features = c('score_ssgw','score_ccat','score_cytotrace')) + theme(legend.position = "right")




COL=c('royalblue3','indianred1','gold1')

pdf('/home/zhangfeng/project/DEV/data/figure/f12_02_score.pdf',width=5,height=7)


############################################################################################
this.obj=readRDS('/home/zhangfeng/project/DEV/data/chickenHeart_GSE149457/D4.rds')
this.xy=GetTissueCoordinates(this.obj)[,c(2,1)]
this.xy[,2]=max(this.xy[,2])-this.xy[,2]

#plot(this.xy)

this_ssgw=this.obj$score_ssgw
this_ccat=this.obj$score_ccat
this_cytotrace=this.obj$score_cytotrace

this_com=.simple_combine(cbind(this_ssgw,this_ccat,this_cytotrace),this.xy)

this_score=this_com$exp_sc_mat1
this_xy=this_com$exp_sc_mat2




##############################
THIS_SCORE=SCORE_SSGW
THIS_COL=.vcol(this_score[,1], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),c(COL))
plot(this_xy,col=THIS_COL,pch=16,main='SSGW')
##############################
THIS_SCORE=SCORE_CCAT
THIS_COL=.vcol(this_score[,2], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),c(COL))
plot(this_xy,col=THIS_COL,pch=16,main='CCAT')
##############################
THIS_SCORE=SCORE_CytoTRACE
THIS_COL=.vcol(this_score[,3], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),c(COL))
plot(this_xy,col=THIS_COL,pch=16,main='CytoTRACE')
#########################################




############################################################################################
this.obj=readRDS('/home/zhangfeng/project/DEV/data/chickenHeart_GSE149457/D7.rds')
this.xy=GetTissueCoordinates(this.obj)[,c(2,1)]
this.xy[,2]=max(this.xy[,2])-this.xy[,2]

#plot(this.xy)

this_ssgw=this.obj$score_ssgw
this_ccat=this.obj$score_ccat
this_cytotrace=this.obj$score_cytotrace

this_com=.simple_combine(cbind(this_ssgw,this_ccat,this_cytotrace),this.xy)

this_score=this_com$exp_sc_mat1
this_xy=this_com$exp_sc_mat2




##############################
THIS_SCORE=SCORE_SSGW
THIS_COL=.vcol(this_score[,1], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),COL)
plot(this_xy,col=THIS_COL,pch=16,main='SSGW')
##############################
THIS_SCORE=SCORE_CCAT
THIS_COL=.vcol(this_score[,2], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),COL)
plot(this_xy,col=THIS_COL,pch=16,main='CCAT')
##############################
THIS_SCORE=SCORE_CytoTRACE
THIS_COL=.vcol(this_score[,3], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),COL)
plot(this_xy,col=THIS_COL,pch=16,main='CytoTRACE')
#########################################




############################################################################################
this.obj=readRDS('/home/zhangfeng/project/DEV/data/chickenHeart_GSE149457/D10.rds')
this.xy=GetTissueCoordinates(this.obj)[,c(2,1)]
this.xy[,2]=max(this.xy[,2])-this.xy[,2]

#plot(this.xy)

this_ssgw=this.obj$score_ssgw
this_ccat=this.obj$score_ccat
this_cytotrace=this.obj$score_cytotrace

this_com=.simple_combine(cbind(this_ssgw,this_ccat,this_cytotrace),this.xy)

this_score=this_com$exp_sc_mat1
this_xy=this_com$exp_sc_mat2




##############################
THIS_SCORE=SCORE_SSGW
THIS_COL=.vcol(this_score[,1], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),COL)
plot(this_xy,col=THIS_COL,pch=16,main='SSGW')
##############################
THIS_SCORE=SCORE_CCAT
THIS_COL=.vcol(this_score[,2], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),COL)
plot(this_xy,col=THIS_COL,pch=16,main='CCAT')
##############################
THIS_SCORE=SCORE_CytoTRACE
THIS_COL=.vcol(this_score[,3], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),COL)
plot(this_xy,col=THIS_COL,pch=16,main='CytoTRACE')
#########################################



############################################################################################
this.obj=readRDS('/home/zhangfeng/project/DEV/data/chickenHeart_GSE149457/D14.rds')
this.xy=GetTissueCoordinates(this.obj)[,c(2,1)]
this.xy[,2]=max(this.xy[,2])-this.xy[,2]

#plot(this.xy)

this_ssgw=this.obj$score_ssgw
this_ccat=this.obj$score_ccat
this_cytotrace=this.obj$score_cytotrace

this_com=.simple_combine(cbind(this_ssgw,this_ccat,this_cytotrace),this.xy)

this_score=this_com$exp_sc_mat1
this_xy=this_com$exp_sc_mat2


##############################
THIS_SCORE=SCORE_SSGW
THIS_COL=.vcol(this_score[,1], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),COL)
plot(this_xy,col=THIS_COL,pch=16,main='SSGW')
##############################
THIS_SCORE=SCORE_CCAT
THIS_COL=.vcol(this_score[,2], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),COL)
plot(this_xy,col=THIS_COL,pch=16,main='CCAT')
##############################
THIS_SCORE=SCORE_CytoTRACE
THIS_COL=.vcol(this_score[,3], c(quantile(THIS_SCORE,0.1),quantile(THIS_SCORE,0.5),quantile(THIS_SCORE,0.9)),COL)
plot(this_xy,col=THIS_COL,pch=16,main='CytoTRACE')
#########################################


plot(rep(1,length(THIS_COL)),col=THIS_COL[order(this_score[,3])],type='h')


dev.off()








