#############################################################################################################
# FitDevo: accurate inference of single-cell developmental potential using sample-specific gene weight
# Author: Feng Zhang
# April 25, 2022
#############################################################################################################


library(Seurat)

.simple_combine <- function(exp_sc_mat1, exp_sc_mat2, FILL=FALSE){
    FILL=FILL
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    ##############################################
    if(FILL==TRUE){
        gene1=rownames(exp_sc_mat)
        gene2=rownames(exp_ref_mat)
        gene12=gene2[which(!gene2 %in% gene1)]
        gene21=gene1[which(!gene1 %in% gene2)]
        exp_sc_mat_add=matrix(0,ncol=ncol(exp_sc_mat),nrow=length(gene12))
        rownames(exp_sc_mat_add)=gene12
        colnames(exp_sc_mat_add)=colnames(exp_sc_mat)
        exp_ref_mat_add=matrix(0,ncol=ncol(exp_ref_mat),nrow=length(gene21))
        rownames(exp_ref_mat_add)=gene21
        colnames(exp_ref_mat_add)=colnames(exp_ref_mat)
        exp_sc_mat=rbind(exp_sc_mat, exp_sc_mat_add)
        exp_ref_mat=rbind(exp_ref_mat, exp_ref_mat_add)
    }
    ############################################
    exp_sc_mat=exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc=rownames(exp_sc_mat)
    gene_ref=rownames(exp_ref_mat)
    gene_over= gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat=exp_sc_mat[which(gene_sc %in% gene_over),]
    exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
    colname_sc=colnames(exp_sc_mat)
    colname_ref=colnames(exp_ref_mat)
    OUT=list()
    OUT$exp_sc_mat1=exp_sc_mat
    OUT$exp_sc_mat2=exp_ref_mat
    OUT$combine=cbind(exp_sc_mat,exp_ref_mat)
    return(OUT)
    }

.toUpper <- function(DATA){
    DATA=DATA
    UPGENE=toupper(rownames(DATA))
    USED_INDEX=which(UPGENE %in% names(which(table(UPGENE)==1)) )
    DATA=DATA[USED_INDEX,]
    rownames(DATA)=UPGENE[USED_INDEX]
    return(DATA)
    }


.normData<-function(mat){
    mat=mat
    mat=mat[which(rownames(mat)!=''),]
    pbmc=CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0, project = "ALL")
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    this_out=as.matrix(pbmc@assays$RNA@data)
    return(this_out)
    }



calScore <- function(MAT, GW){
    MAT=MAT
    GW=GW
    ################################
    COM=.simple_combine(MAT, cbind(GW,GW))
    D1=COM$exp_sc_mat1
    D2=COM$exp_sc_mat2
    D2=D2[,1:(ncol(D2)/2)]
    ####################################
    SCORE=cor(D1, D2)
    return(SCORE)
    }







fitdevo<-function(MAT, BGW, NORM=TRUE, PCNUM=50){
    #################
    library(Seurat)
    #################
    
    #####################################################
    print('FitDevo starts !')
    #####################################################
    print(Sys.time())
    #####################################################
    print('Preprocessing ... ')
    #####################################################

    #####################################################
    # Load data
    MAT=MAT
    BGW=BGW  
    NORM=NORM
    PCNUM=PCNUM

    ##########################################################
    # Solve big matrix
    tooLargeLimit=50000
    tooLargeLimitDelta=10000
    #########################################################
    if(ncol(MAT) > tooLargeLimit){
        set.seed(123)
        splitBy= (seq(ncol(MAT))-1) %/% (tooLargeLimit - tooLargeLimitDelta)
        shuffle_index=shuffle(ncol(MAT))
        splitBy_shuffle=splitBy[shuffle_index] 
        lst = split(colnames(MAT), splitBy_shuffle)
        result_shuffle = unlist(lapply(lst, function(x){fitdevo(MAT[, x], BGW, NORM, PCNUM)}))
        result=result_shuffle
        result[shuffle_index]=result_shuffle
        return(result)    
    }else{
    
    # Preprocess
    this_mat=.toUpper(MAT)
    NNN=min(c(PCNUM,ncol(this_mat)-1))

    #####################################################
    # Create Seurat object
    pbmc=CreateSeuratObject(counts = this_mat, min.cells = 0, min.features = 0, project = "ALL")

    #####################################################
    # Normalization
    if(NORM==TRUE){  
        pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
        }

    #####################################################
    print('Calculating PCs ... ')
    #####################################################
    # Calculate PCs
    pbmc <- FindVariableFeatures(object =pbmc, selection.method = "vst", nfeatures = 2000)
    pbmc <- ScaleData(object = pbmc, features =VariableFeatures(pbmc))
    pbmc <- RunPCA(object = pbmc, npcs=NNN, features = VariableFeatures(pbmc) , ndims.print=1,nfeatures.print=1, seed.use=123)
    PCA=pbmc@reductions$pca@cell.embeddings

    #####################################################
    print('Calculating gene-PC correlation matrix ... ')
    #####################################################
    # Calculate gene-PC correlation matrix
    LOAD=cor(t(as.matrix(pbmc@assays$RNA@data)), PCA) 
    LOAD[which(is.na(LOAD))]=0

    #####################################################
    print('Calculating SSGW ... ')
    #####################################################
    # Use GLM to get pBGW and SSGW
    COM=.simple_combine(cbind(BGW,BGW),LOAD)
    D1=COM$exp_sc_mat1
    D2=COM$exp_sc_mat2
    #################
    Y=D1[,1]
    X=D2
    FIT=glm(Y~.,data=as.data.frame(X),family=binomial(link = "logit"))
    #################
    PRED.Y=predict(FIT,data=as.data.frame(X))
    pBGW=rep(0,length(PRED.Y))
    pBGW[which(PRED.Y>0)]=1
    names(pBGW)=names(PRED.Y)
    SSGW=pBGW+Y
    
    #####################################################
    print('Calculating DP ... ')
    #####################################################
    # Calculate developmental potential (DP) score
    NMAT=as.matrix(pbmc@assays$RNA@data)
    DP=calScore(NMAT, SSGW)
    ###########################

    #######################################
    print('Finished!')
    #######################################
    print(Sys.time())
    #######################################
    return(DP)
    }
    }










