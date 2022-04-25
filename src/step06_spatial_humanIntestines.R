
# /home/toolkit/tools/R4.0.3/bin/R

library(reticulate)
use_python("/home/toolkit/local/bin/python3",required=TRUE)


GW=readRDS(file='/home/zhangfeng/project/DEV/data/RDS/GW.rds')

BGW=GW
BGW[which(GW>0)]=1
BGW[which(GW<=0)]=0
###########################




##########################
setwd('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328')
##########################


library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SeuratDisk)





################################################################
#A1
###############################################################

this.mat=Read10X('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A1/raw_feature_bc_matrix')
this.image=Read10X_Image('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A1/spatial')
DefaultAssay(object = this.image) <- 'Spatial'

this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial')
this.obj[['image']] <- this.image

this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)



plot1 <- VlnPlot(this.obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(this.obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)



##############################################
NDATA=as.matrix(this.obj[['Spatial']]@data)
################################
this.obj$lgr5=NDATA[which(rownames(NDATA)=='LGR5'),]
this.obj$lgr5[which(this.obj$lgr5>0)]=1
#SpatialFeaturePlot(this.obj, features = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1,cols=c('grey90','red1'))
pt1=SpatialPlot(this.obj, group.by = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1.5,cols=c('grey90','red1'))+NoLegend()
pt1


#################################

#################

COUNT=as.matrix(this.obj[['Spatial']]@counts)
RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>0),]
#################
SCORE_SSGW=.runMethod_SSGW(COUNT)
SCORE_CCAT=.runMethod_CCAT(.normData(COUNT))
SCORE_CytoTRACE=.runMethod_CytoTRACE(COUNT)

#########################
this.obj$score_ssgw=SCORE_SSGW
this.obj$score_ccat=SCORE_CCAT
this.obj$score_cytotrace=SCORE_CytoTRACE


#########################################
pt2=SpatialFeaturePlot(this.obj, features = c("score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=1.5,alpha = c(1,1),ncol=3)

pdf('/home/zhangfeng/project/DEV/data/figure/f09_01_A1.pdf',width=12,height=7)
print(pt1)
print(pt2)
dev.off()

#################################
saveRDS(this.obj, file='/home/zhangfeng/project/DEV/data/RDS/Spatial_intestine_A1.rds')
########################







################################################################
#A2
###############################################################

this.mat=Read10X('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A2/raw_feature_bc_matrix')
this.image=Read10X_Image('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A2/spatial')
DefaultAssay(object = this.image) <- 'Spatial'

this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial')
this.obj[['image']] <- this.image

this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)



plot1 <- VlnPlot(this.obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(this.obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)



##############################################
NDATA=as.matrix(this.obj[['Spatial']]@data)
################################
this.obj$lgr5=NDATA[which(rownames(NDATA)=='LGR5'),]
this.obj$lgr5[which(this.obj$lgr5>0)]=1
#SpatialFeaturePlot(this.obj, features = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1,cols=c('grey90','red1'))
pt1=SpatialPlot(this.obj, group.by = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1.5,cols=c('grey90','red1'))+NoLegend()

#################################

#################

COUNT=as.matrix(this.obj[['Spatial']]@counts)
RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>0),]
#################
SCORE_SSGW=.runMethod_SSGW(COUNT)
SCORE_CCAT=.runMethod_CCAT(.normData(COUNT))
SCORE_CytoTRACE=.runMethod_CytoTRACE(COUNT)

#########################
this.obj$score_ssgw=SCORE_SSGW
this.obj$score_ccat=SCORE_CCAT
this.obj$score_cytotrace=SCORE_CytoTRACE


#########################################
pt2=SpatialFeaturePlot(this.obj, features = c("score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=1.5,alpha = c(1,1),ncol=3)

pdf('/home/zhangfeng/project/DEV/data/figure/f09_01_A2.pdf',width=12,height=7)
print(pt1)
print(pt2)
dev.off()

#################################
saveRDS(this.obj, file='/home/zhangfeng/project/DEV/data/RDS/Spatial_intestine_A2.rds')
########################






################################################################
#A3
###############################################################

this.mat=Read10X('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A3/raw_feature_bc_matrix')
this.image=Read10X_Image('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A3/spatial')
DefaultAssay(object = this.image) <- 'Spatial'

this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial')
this.obj[['image']] <- this.image

this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)



plot1 <- VlnPlot(this.obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(this.obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)



##############################################
NDATA=as.matrix(this.obj[['Spatial']]@data)
################################
this.obj$lgr5=NDATA[which(rownames(NDATA)=='LGR5'),]
this.obj$lgr5[which(this.obj$lgr5>0)]=1
#SpatialFeaturePlot(this.obj, features = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1,cols=c('grey90','red1'))
pt1=SpatialPlot(this.obj, group.by = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1.5,cols=c('grey90','red1'))+NoLegend()

#################################

#################

COUNT=as.matrix(this.obj[['Spatial']]@counts)
RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>0),]
#################
SCORE_SSGW=.runMethod_SSGW(COUNT)
SCORE_CCAT=.runMethod_CCAT(.normData(COUNT))
SCORE_CytoTRACE=.runMethod_CytoTRACE(COUNT)

#########################
this.obj$score_ssgw=SCORE_SSGW
this.obj$score_ccat=SCORE_CCAT
this.obj$score_cytotrace=SCORE_CytoTRACE


#########################################
pt2=SpatialFeaturePlot(this.obj, features = c("score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=1.5,alpha = c(1,1),ncol=3)

pdf('/home/zhangfeng/project/DEV/data/figure/f09_01_A3.pdf',width=12,height=7)
print(pt1)
print(pt2)
dev.off()

#################################
saveRDS(this.obj, file='/home/zhangfeng/project/DEV/data/RDS/Spatial_intestine_A3.rds')
########################






################################################################
#A4
###############################################################

this.mat=Read10X('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A4/raw_feature_bc_matrix')
this.image=Read10X_Image('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A4/spatial')
DefaultAssay(object = this.image) <- 'Spatial'

this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial')
this.obj[['image']] <- this.image

this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)



plot1 <- VlnPlot(this.obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(this.obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)



##############################################
NDATA=as.matrix(this.obj[['Spatial']]@data)
################################
this.obj$lgr5=NDATA[which(rownames(NDATA)=='LGR5'),]
this.obj$lgr5[which(this.obj$lgr5>0)]=1
#SpatialFeaturePlot(this.obj, features = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1,cols=c('grey90','red1'))
pt1=SpatialPlot(this.obj, group.by = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1.5,cols=c('grey90','red1'))+NoLegend()

#################################

#################

COUNT=as.matrix(this.obj[['Spatial']]@counts)
RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>0),]
#################
SCORE_SSGW=.runMethod_SSGW(COUNT)
SCORE_CCAT=.runMethod_CCAT(.normData(COUNT))
SCORE_CytoTRACE=.runMethod_CytoTRACE(COUNT)

#########################
this.obj$score_ssgw=SCORE_SSGW
this.obj$score_ccat=SCORE_CCAT
this.obj$score_cytotrace=SCORE_CytoTRACE

#########################################
pt2=SpatialFeaturePlot(this.obj, features = c("score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=1.5,alpha = c(1,1),ncol=3)

pdf('/home/zhangfeng/project/DEV/data/figure/f09_01_A4.pdf',width=12,height=7)
print(pt1)
print(pt2)
dev.off()

#################################
saveRDS(this.obj, file='/home/zhangfeng/project/DEV/data/RDS/Spatial_intestine_A4.rds')
########################





################################################################
#A6
###############################################################

this.mat=Read10X('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A6/raw_feature_bc_matrix')
this.image=Read10X_Image('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A6/spatial')
DefaultAssay(object = this.image) <- 'Spatial'

this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial')
this.obj[['image']] <- this.image

this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)



plot1 <- VlnPlot(this.obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(this.obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)



##############################################
NDATA=as.matrix(this.obj[['Spatial']]@data)
################################
this.obj$lgr5=NDATA[which(rownames(NDATA)=='LGR5'),]
this.obj$lgr5[which(this.obj$lgr5>0)]=1
#SpatialFeaturePlot(this.obj, features = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1,cols=c('grey90','red1'))
pt1=SpatialPlot(this.obj, group.by = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1.5,cols=c('grey90','red1'))+NoLegend()

#################################

#################

COUNT=as.matrix(this.obj[['Spatial']]@counts)
RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>0),]
#################
SCORE_SSGW=.runMethod_SSGW(COUNT)
SCORE_CCAT=.runMethod_CCAT(.normData(COUNT))
SCORE_CytoTRACE=.runMethod_CytoTRACE(COUNT)

#########################
this.obj$score_ssgw=SCORE_SSGW
this.obj$score_ccat=SCORE_CCAT
this.obj$score_cytotrace=SCORE_CytoTRACE


#########################################
pt2=SpatialFeaturePlot(this.obj, features = c("score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=1.5,alpha = c(1,1),ncol=3)

pdf('/home/zhangfeng/project/DEV/data/figure/f09_01_A6.pdf',width=12,height=7)
print(pt1)
print(pt2)
dev.off()

#################################
saveRDS(this.obj, file='/home/zhangfeng/project/DEV/data/RDS/Spatial_intestine_A6.rds')
########################










################################################################
#A7
###############################################################

this.mat=Read10X('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A7/raw_feature_bc_matrix')
this.image=Read10X_Image('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A7/spatial')
DefaultAssay(object = this.image) <- 'Spatial'

this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial')
this.obj[['image']] <- this.image

this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)



plot1 <- VlnPlot(this.obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(this.obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)



##############################################
NDATA=as.matrix(this.obj[['Spatial']]@data)
################################
this.obj$lgr5=NDATA[which(rownames(NDATA)=='LGR5'),]
this.obj$lgr5[which(this.obj$lgr5>0)]=1
#SpatialFeaturePlot(this.obj, features = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1,cols=c('grey90','red1'))
pt1=SpatialPlot(this.obj, group.by = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1.5,cols=c('grey90','red1'))+NoLegend()

#################################

#################

COUNT=as.matrix(this.obj[['Spatial']]@counts)
RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>0),]
#################
SCORE_SSGW=.runMethod_SSGW(COUNT)
SCORE_CCAT=.runMethod_CCAT(.normData(COUNT))
SCORE_CytoTRACE=.runMethod_CytoTRACE(COUNT)

#########################
this.obj$score_ssgw=SCORE_SSGW
this.obj$score_ccat=SCORE_CCAT
this.obj$score_cytotrace=SCORE_CytoTRACE


#########################################
pt2=SpatialFeaturePlot(this.obj, features = c("score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=1.5,alpha = c(1,1),ncol=3)

pdf('/home/zhangfeng/project/DEV/data/figure/f09_01_A7.pdf',width=12,height=7)
print(pt1)
print(pt2)
dev.off()

#################################
saveRDS(this.obj, file='/home/zhangfeng/project/DEV/data/RDS/Spatial_intestine_A7.rds')
########################






################################################################
#A8
###############################################################

this.mat=Read10X('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A8/raw_feature_bc_matrix')
this.image=Read10X_Image('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A8/spatial')
DefaultAssay(object = this.image) <- 'Spatial'

this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial')
this.obj[['image']] <- this.image

this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)



plot1 <- VlnPlot(this.obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(this.obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)



##############################################
NDATA=as.matrix(this.obj[['Spatial']]@data)
################################
this.obj$lgr5=NDATA[which(rownames(NDATA)=='LGR5'),]
this.obj$lgr5[which(this.obj$lgr5>0)]=1
#SpatialFeaturePlot(this.obj, features = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1,cols=c('grey90','red1'))
pt1=SpatialPlot(this.obj, group.by = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1.5,cols=c('grey90','red1'))+NoLegend()

#################################

#################

COUNT=as.matrix(this.obj[['Spatial']]@counts)
RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>0),]
#################
SCORE_SSGW=.runMethod_SSGW(COUNT)
SCORE_CCAT=.runMethod_CCAT(.normData(COUNT))
SCORE_CytoTRACE=.runMethod_CytoTRACE(COUNT)

#########################
this.obj$score_ssgw=SCORE_SSGW
this.obj$score_ccat=SCORE_CCAT
this.obj$score_cytotrace=SCORE_CytoTRACE


#########################################
pt2=SpatialFeaturePlot(this.obj, features = c("score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=1.5,alpha = c(1,1),ncol=3)

pdf('/home/zhangfeng/project/DEV/data/figure/f09_01_A8.pdf',width=12,height=7)
print(pt1)
print(pt2)
dev.off()

#################################
saveRDS(this.obj, file='/home/zhangfeng/project/DEV/data/RDS/Spatial_intestine_A8.rds')
########################







################################################################
#A9
###############################################################

this.mat=Read10X('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A9/raw_feature_bc_matrix')
this.image=Read10X_Image('/home/zhangfeng/project/DEV/data/Human_Intestinal_GSE158328/A9/spatial')
DefaultAssay(object = this.image) <- 'Spatial'

this.obj=CreateSeuratObject(counts =this.mat, project = 'this', assay = 'Spatial')
this.obj[['image']] <- this.image

this.obj <- NormalizeData(object = this.obj, normalization.method = "LogNormalize", scale.factor = 10000)



plot1 <- VlnPlot(this.obj, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(this.obj, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)



##############################################
NDATA=as.matrix(this.obj[['Spatial']]@data)
################################
this.obj$lgr5=NDATA[which(rownames(NDATA)=='LGR5'),]
this.obj$lgr5[which(this.obj$lgr5>0)]=1
#SpatialFeaturePlot(this.obj, features = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1,cols=c('grey90','red1'))
pt1=SpatialPlot(this.obj, group.by = c('lgr5'),pt.size.factor=1.5,alpha = c(1,1),ncol=1.5,cols=c('grey90','red1'))+NoLegend()

#################################

#################

COUNT=as.matrix(this.obj[['Spatial']]@counts)
RSUM=rowSums(COUNT)
COUNT=COUNT[which(RSUM>0),]
#################
SCORE_SSGW=.runMethod_SSGW(COUNT)
SCORE_CCAT=.runMethod_CCAT(.normData(COUNT))
SCORE_CytoTRACE=.runMethod_CytoTRACE(COUNT)

#########################
this.obj$score_ssgw=SCORE_SSGW
this.obj$score_ccat=SCORE_CCAT
this.obj$score_cytotrace=SCORE_CytoTRACE


#########################################
pt2=SpatialFeaturePlot(this.obj, features = c("score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=1.5,alpha = c(1,1),ncol=3)

pdf('/home/zhangfeng/project/DEV/data/figure/f09_01_A9.pdf',width=12,height=7)
print(pt1)
print(pt2)
dev.off()

#################################
saveRDS(this.obj, file='/home/zhangfeng/project/DEV/data/RDS/Spatial_intestine_A9.rds')
########################




cor_NA=function(x,y){ 
    used_index=which((!is.na(x))&(!is.na(y)))
    ccc=cor(x[used_index],y[used_index])
    return(ccc)
    }

MAT=matrix(0,ncol=8,nrow=3)

this.obj=A1
this.obj$lgr5=this.obj[['Spatial']]@data[which(rownames(this.obj[['Spatial']]@data)=='LGR5'),]



MAT[1,1]=cor_NA(this.obj$lgr5,this.obj$score_ssgw)
MAT[2,1]=cor_NA(this.obj$lgr5,this.obj$score_ccat)
MAT[3,1]=cor_NA(this.obj$lgr5,this.obj$score_cytotrace)


this.obj=A2





this.obj$lgr5=this.obj[['Spatial']]@data[which(rownames(this.obj[['Spatial']]@data)=='LGR5'),]
MAT[1,2]=cor_NA(this.obj$lgr5,this.obj$score_ssgw)
MAT[2,2]=cor_NA(this.obj$lgr5,this.obj$score_ccat)
MAT[3,2]=cor_NA(this.obj$lgr5,this.obj$score_cytotrace)

this.obj=A3
this.obj$lgr5=this.obj[['Spatial']]@data[which(rownames(this.obj[['Spatial']]@data)=='LGR5'),]

MAT[1,3]=cor_NA(this.obj$lgr5,this.obj$score_ssgw)
MAT[2,3]=cor_NA(this.obj$lgr5,this.obj$score_ccat)
MAT[3,3]=cor_NA(this.obj$lgr5,this.obj$score_cytotrace)

this.obj=A4
this.obj$lgr5=this.obj[['Spatial']]@data[which(rownames(this.obj[['Spatial']]@data)=='LGR5'),]

MAT[1,4]=cor_NA(this.obj$lgr5,this.obj$score_ssgw)
MAT[2,4]=cor_NA(this.obj$lgr5,this.obj$score_ccat)
MAT[3,4]=cor_NA(this.obj$lgr5,this.obj$score_cytotrace)

this.obj=A6
this.obj$lgr5=this.obj[['Spatial']]@data[which(rownames(this.obj[['Spatial']]@data)=='LGR5'),]

MAT[1,5]=cor_NA(this.obj$lgr5,this.obj$score_ssgw)
MAT[2,5]=cor_NA(this.obj$lgr5,this.obj$score_ccat)
MAT[3,5]=cor_NA(this.obj$lgr5,this.obj$score_cytotrace)

this.obj=A7
this.obj$lgr5=this.obj[['Spatial']]@data[which(rownames(this.obj[['Spatial']]@data)=='LGR5'),]

MAT[1,6]=cor_NA(this.obj$lgr5,this.obj$score_ssgw)
MAT[2,6]=cor_NA(this.obj$lgr5,this.obj$score_ccat)
MAT[3,6]=cor_NA(this.obj$lgr5,this.obj$score_cytotrace)

this.obj=A8
this.obj$lgr5=this.obj[['Spatial']]@data[which(rownames(this.obj[['Spatial']]@data)=='LGR5'),]

MAT[1,7]=cor_NA(this.obj$lgr5,this.obj$score_ssgw)
MAT[2,7]=cor_NA(this.obj$lgr5,this.obj$score_ccat)
MAT[3,7]=cor_NA(this.obj$lgr5,this.obj$score_cytotrace)

this.obj=A9
this.obj$lgr5=this.obj[['Spatial']]@data[which(rownames(this.obj[['Spatial']]@data)=='LGR5'),]

MAT[1,8]=cor_NA(this.obj$lgr5,this.obj$score_ssgw)
MAT[2,8]=cor_NA(this.obj$lgr5,this.obj$score_ccat)
MAT[3,8]=cor_NA(this.obj$lgr5,this.obj$score_cytotrace)

MAT=t(MAT)
colnames(MAT)=c('SSGW','CCAT','CytoTRACE')
rownames(MAT)=paste0('S',c(1:8))



library('ComplexHeatmap')
library('circlize')


mat=MAT-apply(MAT,1,min)

Heatmap3D(mat)


pdf('/home/zhangfeng/project/DEV/data/figure/f09_02_corHeat_LGR5.pdf',width=3,height=3)

color_fun_3 =colorRamp2(c(0,0.01,0.04 ), c('grey95','indianred1','gold1'))
Heatmap(mat,row_title='',name="d",
        cluster_columns=FALSE, cluster_rows=FALSE,
	      show_column_dend = FALSE, show_row_dend = FALSE, 
	      show_column_names=FALSE, show_row_names=FALSE,
	      col=color_fun_3, border = TRUE  ,
          cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.3f", MAT[i, j]), x, y, gp = gpar(fontsize = 10))
          }  
	      )
dev.off()





pdf('/home/zhangfeng/project/DEV/data/figure/f09_03_SFeatures.pdf',width=16,height=7)
SpatialFeaturePlot(A1, features = c('LGR5',"score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=2,alpha = c(0.1,1),ncol=4,stroke=0)
SpatialFeaturePlot(A2, features = c('LGR5',"score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=2,alpha = c(0.1,1),ncol=4,stroke=0)
SpatialFeaturePlot(A3, features = c('LGR5',"score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=2,alpha = c(0.1,1),ncol=4,stroke=0)
SpatialFeaturePlot(A4, features = c('LGR5',"score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=2,alpha = c(0.1,1),ncol=4,stroke=0)
SpatialFeaturePlot(A6, features = c('LGR5',"score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=2,alpha = c(0.1,1),ncol=4,stroke=0)
SpatialFeaturePlot(A7, features = c('LGR5',"score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=2,alpha = c(0.1,1),ncol=4,stroke=0)
SpatialFeaturePlot(A8, features = c('LGR5',"score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=2,alpha = c(0.1,1),ncol=4,stroke=0)
SpatialFeaturePlot(A9, features = c('LGR5',"score_ssgw","score_ccat","score_cytotrace"),pt.size.factor=2,alpha = c(0.1,1),ncol=4,stroke=0)
dev.off()



