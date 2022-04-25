

# /home/toolkit/tools/R4.0.3/bin/R

source('https://gitee.com/jumphone/public/raw/master/fitdevo.R')
GW=readRDS(url('https://gitee.com/jumphone/public/raw/master/GW.rds'))
BGW=readRDS(url('https://gitee.com/jumphone/public/raw/master/BGW.rds'))


##########################
setwd('/home/zhangfeng/project/DEV/data/CCAT')
##########################


.entriz2symbol <- function(entriz){
    entriz=entriz
    library('org.Hs.eg.db')
    symbol=mapIds(org.Hs.eg.db, entriz, 'SYMBOL','ENTREZID')
    return(symbol)
    }
    
   

.evalMethod_28data <-function(.runMethod, NORM=TRUE){
    library(pROC)
    COR=c()
    AUC=c()
    i=1
    while(i<=28){
        ID=i
        this_file=paste0('/home/zhangfeng/project/DEV/data/CCAT/',ID,'.rds')
        this_data=readRDS(this_file)
        this_mat=this_data$mat
        if(NORM==TRUE){
            this_mat= .normData(this_mat)
            }
        this_tag= this_data$tag
        this_mat=as.matrix(this_mat)
        ####################################
        CSUM=rowSums(t(this_mat))
        RSUM=rowSums(this_mat)
        this_mat=this_mat[which(RSUM>0),]
        ##################################### 
        used_index=which(!is.na(this_tag) & CSUM>0)
        this_mat=this_mat[,used_index]
        this_tag=this_tag[used_index]
        ############################################################################
        #####################################
        this_mat=as.matrix(this_mat)
        #####################################
        this_mat=.toUpper(this_mat)
        this_score=.runMethod(this_mat)
        ###########################
        ##########################
        used_index=which((!is.na(this_score)) )
        ##########################
        this_pred=this_score[used_index]
        this_resp=this_tag[used_index]
        this_roc=roc(response=this_resp, predictor=this_pred, quiet =TRUE,direction='<')
        this_auc=auc(this_roc)[1]
        #######################
        this_cor=cor(this_score[used_index], this_tag[used_index],method='spearman')
        ####################
        print(this_auc)
        print(this_cor)
        ##################
        COR=c(COR, this_cor)
        AUC=c(AUC, this_auc)
        ########################
        print(i)
        i=i+1
        }
    ##############
    OUT=list()
    OUT$cor=COR
    OUT$auc=AUC
    return(OUT)
    }



##################################################
# GW
###################################################
.runMethod_GW<-function(mat){
        this_score=calScore(mat, GW)
        return(this_score)
        }

EVAL_GW_28data=.evalMethod_28data(.runMethod_GW,NORM=TRUE)


##################################################
# BGW
###################################################
.runMethod_BGW<-function(mat){
        this_score=calScore(mat, BGW)
        return(this_score)
        }

EVAL_BGW_28data=.evalMethod_28data(.runMethod_BGW,NORM=TRUE)


##################################################
# FitDevo (SSGW)
###################################################
.runMethod_SSGW <-function(mat){ # FitDevo
        this_score=fitdevo(MAT=mat, BGW=BGW, NORM=FALSE, PCNUM=50 ) 
        return(this_score)
        }

EVAL_SSGW_28data=.evalMethod_28data(.runMethod_SSGW,NORM=TRUE)


##################################################
# CytoTRACE
###################################################
source('/home/zhangfeng/project/DEV/data/CytoTRACE/CytoTRACE.R')

.runMethod_CytoTRACE<-function(mat){
    this_score=CytoTRACE(mat)$CytoTRACE
    return(this_score)
    }

EVAL_CytoTRACE_28data=.evalMethod_28data(.runMethod_CytoTRACE, NORM=FALSE) # CytoTRACE has normalization process


############################################################
#CCAT
#############################################################

library(SCENT)
library('org.Hs.eg.db')
########################
data(net17Jan16)
NET= net17Jan16.m
NET.SYM=.entriz2symbol(rownames(NET))
USED_ROW=which(!is.na(NET.SYM))
NET=NET[USED_ROW,]
rownames(NET)=NET.SYM[USED_ROW]
NET=.toUpper(NET)
##############

CompCCAT <- function (exp.m, ppiA.m){
    if (max(exp.m) > 100) {
        exp.m <- log2(exp.m + 1)
    }
    classMATRIX <- class(exp.m)
    commonEID.v <- intersect(rownames(ppiA.m), rownames(exp.m))
    k.v <- rowSums(ppiA.m[match(commonEID.v, rownames(ppiA.m)),
        ])
    if (classMATRIX == "matrix") {
        ccat.v <- as.vector(cor(exp.m[match(commonEID.v, rownames(exp.m)),
            ], k.v))
    }
    else if (classMATRIX == "dgCMatrix") {
        ccat.v <- as.vector(corSparse(exp.m[match(commonEID.v,
            rownames(exp.m)), ], Matrix(matrix(k.v, ncol = 1))))
    }
    return(ccat.v)
}

.runMethod_CCAT<-function(mat){
    mat=.toUpper(mat)
    this_score=CompCCAT(mat, NET)
    return(this_score)
    }

EVAL_CCAT_28data=.evalMethod_28data(.runMethod_CCAT,NORM=TRUE)





####################################
# Spearman Correlation Coefficient (SCC) 
###################################


###############
# All
###############

print(mean(EVAL_GW_28data$cor))
# 0.6631595
print(mean(EVAL_BGW_28data$cor))
# 0.6769669
print(mean(EVAL_SSGW_28data$cor))
# 0.6942443
print(mean(EVAL_CytoTRACE_28data$cor))
# 0.5940493
print(mean(EVAL_CCAT_28data$cor))
# 0.6256122

###############
# Overlapped set
###############
OVERLAP_INDEX=c(5,6,7,8,13,16,19,21,22,23)

print(mean(EVAL_GW_28data$cor[OVERLAP_INDEX]))
# 0.7229761
print(mean(EVAL_BGW_28data$cor[OVERLAP_INDEX]))
# 0.7410273
print(mean(EVAL_SSGW_28data$cor[OVERLAP_INDEX]))
# 0.7532628
print(mean(EVAL_CytoTRACE_28data$cor[OVERLAP_INDEX]))
# 0.6070611
print(mean(EVAL_CCAT_28data$cor[OVERLAP_INDEX]))
# 0.6662716

###############
# Novel set
###############
NOVEL_INDEX=c(1,2,3,4,9,10,11,12,14,15,17,18,20,24,25,26,27,28)

print(mean(EVAL_GW_28data$cor[NOVEL_INDEX]))
# 0.629928
print(mean(EVAL_BGW_28data$cor[NOVEL_INDEX]))
# 0.6413778
print(mean(EVAL_SSGW_28data$cor[NOVEL_INDEX]))
# 0.6614562
print(mean(EVAL_CytoTRACE_28data$cor[NOVEL_INDEX]))
# 0.5868205
print(mean(EVAL_CCAT_28data$cor[NOVEL_INDEX]))
# 0.6030236






####################################
# Area under the curve (AUC)
###################################

###############
# All
###############
print(mean(EVAL_GW_28data$auc))
# 0.911352
print(mean(EVAL_BGW_28data$auc))
# 0.9198596
print(mean(EVAL_SSGW_28data$auc))
# 0.930488
print(mean(EVAL_CytoTRACE_28data$auc))
# 0.8677143
print(mean(EVAL_CCAT_28data$auc))
# 0.8890971

###############
# Overlapped set
###############
OVERLAP_INDEX=c(5,6,7,8,13,16,19,21,22,23)

print(mean(EVAL_GW_28data$auc[OVERLAP_INDEX]))
# 0.9241374
print(mean(EVAL_BGW_28data$auc[OVERLAP_INDEX]))
# 0.9349498
print(mean(EVAL_SSGW_28data$auc[OVERLAP_INDEX]))
# 0.9422285
print(mean(EVAL_CytoTRACE_28data$auc[OVERLAP_INDEX]))
# 0.8527953
print(mean(EVAL_CCAT_28data$auc[OVERLAP_INDEX]))
# 0.8907045

###############
# Novel set
###############
NOVEL_INDEX=c(1,2,3,4,9,10,11,12,14,15,17,18,20,24,25,26,27,28)

print(mean(EVAL_GW_28data$auc[NOVEL_INDEX]))
# 0.9042491
print(mean(EVAL_BGW_28data$auc[NOVEL_INDEX]))
# 0.9114762
print(mean(EVAL_SSGW_28data$auc[NOVEL_INDEX]))
# 0.9239654
print(mean(EVAL_CytoTRACE_28data$auc[NOVEL_INDEX]))
# 0.8760026
print(mean(EVAL_CCAT_28data$auc[NOVEL_INDEX]))
# 0.8882041






