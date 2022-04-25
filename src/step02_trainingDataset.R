

# /home/toolkit/tools/R4.0.3/bin/R

source('https://gitee.com/jumphone/public/raw/master/fitdevo.R')
GW=readRDS(url('https://gitee.com/jumphone/public/raw/master/GW.rds'))
BGW=readRDS(url('https://gitee.com/jumphone/public/raw/master/BGW.rds'))



##########################
setwd('/home/zhangfeng/project/DEV/data/CytoTRACE')
##########################


.entriz2symbol <- function(entriz){
    entriz=entriz
    library('org.Hs.eg.db')
    symbol=mapIds(org.Hs.eg.db, entriz, 'SYMBOL','ENTREZID')
    return(symbol)
    }
    
    

.evalMethod_17data <-function(.runMethod, NORM=TRUE){

    COR=c()
    i=1
    while(i<=17){
        ID=i
        this_file=paste0('/home/zhangfeng/project/DEV/data/CytoTRACE/',ID,'.rds')
        this_data=readRDS(this_file)
        this_mat=this_data$mat
        if(NORM==TRUE){
            this_mat= .normData(this_mat)
            }
        this_tag= this_data$tag
        #################################### 
        used_index=which(!is.na(this_tag))
        this_mat=this_mat[,used_index]
        this_tag=this_tag[used_index]
        ############################################################################
        #####################################
        this_mat=as.matrix(this_mat)
        #####################################
        this_mat=.toUpper(this_mat)
        this_mat=this_mat[which(rownames(this_mat)!=''),]
        this_score=.runMethod(this_mat)
        ##########################
        used_index=which((!is.na(this_score)) )
        ##########################
        this_pred=this_score[used_index]
        this_resp=this_tag[used_index]
        #######################
        this_cor=cor(this_pred, this_resp,method='spearman')
        COR=c(COR, this_cor)
        ########################
        print(i)
        i=i+1
        }
    ##############
    OUT=list()
    OUT$cor=COR
    return(OUT)
    }

##################################################
# GW
###################################################
.runMethod_GW<-function(mat){
        this_score=calScore(mat, GW)
        return(this_score)
        }

EVAL_GW_17data=.evalMethod_17data(.runMethod_GW,NORM=TRUE)
print(mean(EVAL_GW_17data$cor))
# 0.6347399

##################################################
# BGW
###################################################
.runMethod_BGW<-function(mat){
        this_score=calScore(mat, BGW)
        return(this_score)
        }

EVAL_BGW_17data=.evalMethod_17data(.runMethod_BGW,NORM=TRUE)
print(mean(EVAL_BGW_17data$cor))
# 0.6470112

##################################################
# FitDevo (SSGW)
###################################################
.runMethod_SSGW <-function(mat){ # FitDevo
        this_score=fitdevo(MAT=mat, BGW=BGW, NORM=FALSE, PCNUM=50 ) 
        return(this_score)
        }

EVAL_SSGW_17data=.evalMethod_17data(.runMethod_SSGW,NORM=TRUE)
print(mean(EVAL_SSGW_17data$cor))
# 0.6694486

##################################################
# CytoTRACE
###################################################
source('/home/zhangfeng/project/DEV/data/CytoTRACE/CytoTRACE.R')

.runMethod_CytoTRACE<-function(mat){
    this_score=CytoTRACE(mat)$CytoTRACE
    return(this_score)
    }

EVAL_CytoTRACE_17data=.evalMethod_17data(.runMethod_CytoTRACE, NORM=FALSE) # CytoTRACE has normalization process
print(mean(EVAL_CytoTRACE_17data$cor))
# 0.6046702




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

EVAL_CCAT_17data=.evalMethod_17data(.runMethod_CCAT,NORM=TRUE)
print(mean(EVAL_CCAT_17data$cor))
# 0.5011177  





















