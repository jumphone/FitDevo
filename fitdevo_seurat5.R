#############################################################################################################
# FitDevo: accurate inference of single-cell developmental potential using sample-specific gene weight
# Author: Feng Zhang
# Date of first version (v1.0.0): April 25, 2022
#############################################################################################################
print('#############################################################')
print('')
print('Thanks for using FitDevo v1.2')
print('')
print('#############################################################')
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
    #this_out=as.matrix(pbmc@assays$RNA@data)
    this_out=pbmc@assays$RNA@data
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





fitdevo<-function(MAT, BGW, NORM=TRUE, PCNUM=50, VARGENE=2000, tooLargeLimit=50000, SEED=123, SS=TRUE, DETAIL=FALSE){
    #################
    library(Seurat)
    library(qlcMatrix)
    #################

    #####################################################
    # Load data
    MAT=MAT
    BGW=BGW  
    NORM=NORM
    PCNUM=PCNUM
    VARGENE=VARGENE
    SS=SS
    DETAIL=DETAIL
    SEED=SEED
    ##########################################################
    # Solve big matrix
    tooLargeLimit=tooLargeLimit
        
    #####################################################
    print('FitDevo starts !')
    #####################################################
    print(Sys.time())
    #####################################################
    print('Preprocessing ... ')
    #####################################################
    
    # Preprocess
    this_mat=.toUpper(MAT)

    #############################################
    # Create Seurat object
    pbmc=CreateSeuratObject(counts = this_mat, min.cells = 0, min.features = 0, project = "ALL")

    #############################################
    # Normalization
    if(NORM==TRUE){  
        pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000) 
        }

    ##############################################   
    NMAT=pbmc@assays$RNA@data
    NCOL=ncol(NMAT)
    USED_GW=BGW
    USED_GW=USED_GW[which(names(USED_GW) %in% rownames(NMAT))]
    
    ################################
    if(length(unique(USED_GW))==1){stop('GW is not suitable !!!')}
    ################################

    if(SS==TRUE){

        ###################################################################
        ###################################################################
        # Calculate SSGW - START

        ################################################
        if(NCOL > tooLargeLimit){
            ############################
            print('Cell number is larger than tooLargeLimit. Conduct down-sampling...')
            print('Cell number:')
            print(NCOL)
            print('tooLargeLimit:')
            print(tooLargeLimit)
            ############################
            set.seed(SEED)
            used_index=sample(c(1:NCOL), tooLargeLimit)
            pbmc=subset(pbmc, cells=colnames(pbmc)[used_index])
            }

        
        #####################################################
        print('Calculating PCs ... ')
        #####################################################
        # Calculate PCs
        pbmc <- FindVariableFeatures(object =pbmc, selection.method = "vst", nfeatures = VARGENE)
        pbmc <- ScaleData(object = pbmc, features =VariableFeatures(pbmc))
        ####################################
        NNN=min(c(PCNUM,ncol(pbmc)-1))
        pbmc <- RunPCA(object = pbmc, npcs=NNN, features = VariableFeatures(pbmc) , ndims.print=1,nfeatures.print=1, seed.use=SEED)
        ####################################
        PCA=pbmc@reductions$pca@cell.embeddings

        #####################################################
        print('Calculating gene-PC correlation matrix ... ')
        #####################################################
        # Calculate gene-PC correlation matrix
        options(warn=-1) 
        LOAD=Matrix(corSparse(Matrix::t(pbmc@assays$RNA@data), Matrix(PCA)))
        rownames(LOAD)=rownames(pbmc@assays$RNA@data)
        colnames(LOAD)=colnames(PCA)
        options(warn=1)
        LOAD[which(is.na(LOAD))]=0

        #####################################################
        print('Calculating SSGW ... ')
        #####################################################
        # Use GLM to get pBGW and SSGW
        COM=.simple_combine(cbind(USED_GW, USED_GW),LOAD)
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
        USED_GW=SSGW
        #######################
        # Calculate SSGW - END
        ###################################################################
        ###################################################################

        }

    #####################################################
    print('Calculating DP ... ')
    #####################################################
    # Calculate developmental potential (DP) score
    options(warn=-1)
    DP=as.vector(corSparse(NMAT[match(names(USED_GW), rownames(NMAT)), ], Matrix(matrix(USED_GW, ncol = 1))))
    options(warn=1)
    ###########################

    #######################################
    print('Finished!')
    #######################################
    print(Sys.time())
    #######################################
    if(DETAIL==TRUE){

        OUT=list()
        OUT$DP=DP
        OUT$USED_GW=USED_GW
        OUT$parameters=list()
        OUT$parameters$NORM=NORM
        OUT$parameters$PCNUM=PCNUM
        OUT$parameters$VARGENE=VARGENE
        OUT$parameters$SS=SS
        OUT$parameters$DETAIL=DETAIL
        OUT$parameters$tooLargeLimit=tooLargeLimit
        OUT$parameters$SEED=SEED
        return(OUT)

        }else{

        return(DP)

        }
    ##########################################
    }















###############################################################################################################################

# 2022.07.01

###############################################################################################################################

.norm1<-function(x){
    if(min(x)!=max(x)){
        y=(x-min(x))/(max(x)-min(x))
       }else{
        y=x   
       }
    return(y)
    }


#########################################################################################


.vcol<-function(VALUE, CV, CN){
    VALUE=.norm1(VALUE)
    CV=.norm1(CV)
    CN=CN
    if(max(VALUE)==min(VALUE)){return(rep('grey70',length(VALUE)))}
    ##############################################
    CMAT=c()
    i=1
    while(i<=length(CN)){
        CMAT=cbind(CMAT, col2rgb(CN[i]))
        i=i+1
        }
     CMAT=t(CMAT)
     #################################################
     RED.V=CMAT[,1]
     GRE.V=CMAT[,2]
     BLU.V=CMAT[,3]
     ###############################################
     options(warn=-1)
     #####################
     train.data=data.frame(x=CV)
     fit.red=loess(RED.V~x,data=train.data)
     fit.gre=loess(GRE.V~x,data=train.data)
     fit.blu=loess(BLU.V~x,data=train.data)
     ###############################################
     pred.data=data.frame(x=VALUE)
     pred.red=predict(fit.red,newdata=pred.data)
     pred.gre=predict(fit.gre,newdata=pred.data)
     pred.blu=predict(fit.blu,newdata=pred.data)
     ####################
     options(warn=1)
     ##############################################
     .filterCol<-function(x){
         y=x
         y[which(x<0)]=0
         y[which(x>255)]=255
         return(y)
         }
     ####################
     pred.red=.filterCol(pred.red)
     pred.gre=.filterCol(pred.gre)
     pred.blu=.filterCol(pred.blu)

     ########################
     PRED.COL=rgb(red=pred.red/255,green=pred.gre/255,blue=pred.blu/255)
     return(PRED.COL)
    }



#########################################################################################


############################################################
# "cart2clock" and "clock2cart" refer to "IDPmisc" package
############################################################

.getXY = function (x, y = NULL, unidim.allowed = TRUE){
    if (missing(x))
        xarg <- "x"
    else xarg <- deparse(substitute(x))
    if (missing(y))
        yarg <- "y"
    else yarg <- deparse(substitute(y))
    if (is.matrix(x) | is.data.frame(x)) {
        if (ncol(x) > 1) {
            if (is.null(y))
                y <- x[, 2]
            else stop("'", xarg, "' must have only 1 column when y is supplied                                                                                                                                                               separately\n")
        }
        x <- x[, 1]
    }
    else if (is.list(x)) {
        if (length(x) > 1) {
            y <- x[[2]]
            x <- x[[1]]
            if (length(y) != length(x))
                stop("First and second element of the list must have identical                                                                                                                                                               length!\n")
        }
        else x <- x[[1]]
    }
    if (is.null(y)) {
        if (!unidim.allowed)
            stop("'", yarg, "' is not defined!\n")
        y <- x
        x <- 1:length(x)
    }
    else {
        y <- unlist(y)
        if (length(y) != length(x))
            stop("Vector '", yarg, "' and 'x' must have identical lengths!\n")
    }
    return(data.frame(x = x, y = y))
}



.cart2clock = function (x, y = NULL, circle=360){
    xy <- .getXY(x, y, unidim.allowed = FALSE)
    x <- xy$x
    y <- xy$y
    phi <- (atan(x/y)/2/pi * circle + ifelse(y >= 0, circle,
        1.5 * circle))%%circle
    return(data.frame(rho = sqrt(x * x + y * y), phi = phi))
    }

 
.clock2cart = function (rho, phi = NULL, circle=360){
    xy <- .getXY(rho, phi, unidim.allowed = FALSE)
    rho <- xy$x
    phi <- xy$y
    return(data.frame(x = rho * sin(phi/circle * 2 * pi), y = rho *
        cos(phi/circle * 2 * pi)))
    }
##################################################################################################


.buildGrid <- function(VEC, N=30, SHOW=FALSE, COL='grey70'){
    ###############################
    VEC.E=VEC
    COL=COL
    SHOW=SHOW
    delta=0.000001
    N=N
    #################################
    if(SHOW==TRUE){
        plot(VEC.E,col=COL,pch=16,cex=0.2)
        }
    #############################
    X=VEC[,1]
    Y=VEC[,2]

    ########################################################
    this_step_x=( max(X) - min(X) - delta)/N
    this_step_y=( max(Y) - min(Y) - delta)/N
    NUM_CUT=1
    #######################################################
    INDEX_LIST=list()
    CENTER_LIST=list()
    #######################################################
    this_x=min(X)
    while(this_x<max(X)){
        this_y=min(Y)
        while(this_y<max(Y)){  
            #######################################################
            this_in_index = which(VEC.E[,1]>=this_x &  VEC.E[,1] <this_x+this_step_x &
                                  VEC.E[,2]>=this_y &  VEC.E[,2] <this_y+this_step_y)                                         
            this_center=c(this_x+this_step_x/2,this_y+this_step_y/2)                  
            if(length(this_in_index)>=NUM_CUT){
                #######################################################
                INDEX_LIST=c(INDEX_LIST, list(this_in_index))
                CENTER_LIST=c(CENTER_LIST, list(this_center))
                #######################################################
                if(SHOW==TRUE){
                    points(this_center[1],this_center[2],col= 'black',pch=16,cex=0.5)
                    }
                ######################
                }
                 
            this_y=this_y+this_step_y
            }
        this_x=this_x+this_step_x
        }
    ############################################
    OUT=list()
    OUT$VEC=VEC.E
    OUT$INDEX_LIST=INDEX_LIST
    OUT$CENTER_LIST=CENTER_LIST
    OUT$this_step_x=this_step_x
    OUT$this_step_y=this_step_y
    ########################################
    return(OUT)
    }





#########################################################################################


.buildNet<-function(OUT, CUT=1, SHOW=FALSE, COL='grey70'){
    #########################################
    #library(igraph)
    library(stringr)
    ###########################################
    OUT=OUT
    SHOW=SHOW
    CUT=CUT
    COL=COL
    VEC=OUT$VEC
    INDEX_LIST=OUT$INDEX_LIST
    CENTER_LIST=OUT$CENTER_LIST
    this_step_x=OUT$this_step_x
    this_step_y=OUT$this_step_y
    delta=0.00001
    ###############################
    if(SHOW==TRUE){
        plot(VEC,col=COL,pch=16,cex=0.2)
        }
    #############################
    CENTER_VEC=c()
    i=1
    while(i<=length(CENTER_LIST)){
        CENTER_VEC=cbind(CENTER_VEC, CENTER_LIST[[i]])
        i=i+1}
    CENTER_VEC=t(CENTER_VEC)
    OUT$CENTER_VEC=CENTER_VEC
    ##############################
    p1=c()
    p2=c()
    
    #############################
    CNUM=length(CENTER_LIST)
    i=1
    while(i<=CNUM){
    
        this_p1_loc=CENTER_LIST[[i]] 
        this_p1=paste0('P',as.character(i))
    
        used_j=which( ( abs(CENTER_VEC[,1]-this_p1_loc[1]) <= this_step_x+delta) 
                     & ( abs(CENTER_VEC[,2]-this_p1_loc[2]) <= this_step_y+delta) )
        
        for(j in used_j){
            this_p2=paste0('P',as.character(j))
            
            this_p2_loc=CENTER_LIST[[j]]
 
            ######################
            
            if(length(INDEX_LIST[[i]])>=CUT & 
               length(INDEX_LIST[[j]])>=CUT &
               this_p1 != this_p2){
                
                p1=c(p1,this_p1)
                p2=c(p2,this_p2) 
                if(SHOW==TRUE){
                    segments(x0=this_p1_loc[1],x1=this_p2_loc[1],
                           y0=this_p1_loc[2],y1=this_p2_loc[2],
                       
                           col='black')
                    } 
                }
    
            ############
            }
       
            #if(i%%10==1){print(paste0(i,'/',CNUM))}
            i=i+1
       }
    ##########################  
    
    
    
    TAG=c()
    i=1
    while(i<=length(p1)){
        this_p1=p1[i]
        this_p2=p2[i]
        sorted_pair=sort(c(this_p1,this_p2))
        this_tag=paste0(sorted_pair[1],'|',sorted_pair[2])
        TAG=c(TAG, this_tag)
        i=i+1}
    TAG=unique(TAG)
    
    p1=c()
    p2=c()
    i=1
    while(i<=length(TAG)){
        this_p1=strsplit(TAG[i],'\\|')[[1]][1]
        this_p2=strsplit(TAG[i],'\\|')[[1]][2]
        p1=c(p1,this_p1)
        p2=c(p2,this_p2)
        i=i+1}
    
    ##########################
    
    #library(igraph)
    OUT$p1=p1
    OUT$p2=p2
    NET = cbind(p1,p2) 
    g <- igraph::make_graph(t(NET),directed = FALSE)
    ALLNAME=paste0('P',1:CNUM)
    ADD=ALLNAME[which(! ALLNAME %in% igraph::as_ids(igraph::V(g)))]
    g <- g + ADD
    
    ##########################
    DIST=igraph::distances(g, v = igraph::V(g), to = igraph::V(g), mode = c("all"))
    library(stringr)
    DIST.NUM=as.numeric(str_replace(colnames(DIST),'P',''))
    DIST=DIST[order(DIST.NUM),order(DIST.NUM)]
    ###########################
    #library(igraph)
    CPT=igraph::components(g)
    MAXC=which(CPT$csize==max(CPT$csize))[1]
    ############
    library(stringr)
    ###############
    USED_NAME=names(which(CPT$membership==MAXC))
    USED=as.numeric(str_replace(USED_NAME,'P',''))
    
    USED_NAME=USED_NAME[order(USED)]
    USED=USED[order(USED)]
    if(SHOW==TRUE){
        points(OUT$CENTER_VEC[USED,],col='red',pch=16, cex=0.5)
        }
    USED_INDEX=c()
    i=1
    while(i<=length(USED)){
        USED_INDEX=c(USED_INDEX,INDEX_LIST[[USED[i]]])
        i=i+1}
    
    
    ######################################
    POINT_INDEX=rep(NA, nrow(VEC))
    i=1
    while(i<=length(INDEX_LIST)){
        POINT_INDEX[INDEX_LIST[[i]]]=i
        i=i+1
        }


    ###################
    OUT$POINT.INDEX=POINT_INDEX
    ###################
    OUT$GRAPH=g
    OUT$DIST=DIST
    OUT$USED=USED
    OUT$USED_NAME=USED_NAME
    OUT$USED_INDEX=USED_INDEX
    return(OUT)
    ###########################
    }




.gridValue <- function(OUT, SHOW=FALSE){
    OUT=OUT
    INDEX_LIST=OUT$INDEX_LIST
    VALUE=OUT$VALUE
    SHOW=SHOW
    USED=OUT$USED
    ######################
    
    ################
    CENTER_VALUE=c()
    i=1
    while(i<=length(INDEX_LIST)){
        this_value = mean(VALUE[INDEX_LIST[[i]]])
        CENTER_VALUE=c(CENTER_VALUE,this_value)
        i=i+1}
    ############
    CENTER_VEC=OUT$CENTER_VEC

    ##################
    VALUE.COL=.vcol(.norm1(VALUE), c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    ####################
    
    #################################
    CENTER.GRID.COL=.vcol(.norm1(CENTER_VALUE), c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))

    ################
    if(SHOW==TRUE){   
        plot(OUT$VEC,col='grey80',pch=16,cex=0.5)
        points(OUT$CENTER_VEC[USED,],col=CENTER.GRID.COL[USED],pch=15, cex=1.5)
        }

    ######################
    CENTER.POINT.COL=rep('grey70',nrow(OUT$VEC))
    i=1
    while(i<=length(USED)){
        this_index=INDEX_LIST[[USED[i]]]
        CENTER.POINT.COL[this_index]=CENTER.GRID.COL[USED][i]
        i=i+1}

    ################
    OUT$POINT.COL=VALUE.COL
    OUT$POINT.CENTER.COL=CENTER.POINT.COL    

    OUT$CENTER.VALUE=CENTER_VALUE     
    OUT$CENTER.GRID.COL=CENTER.GRID.COL
    ##################
    return(OUT)
    }




.drawArrow <- function(OUT, P=0.9, SHOW=TRUE, COL='grey70', CEX=0.5, LWD=1.5, BD=TRUE, AC='grey10',AL=0.4){
    
    ################
    CEX=CEX
    LWD=LWD
    BD=BD
    AC=AC
    AL=AL
    OUT=OUT
    SHOW=SHOW
    P=P
    COL=COL
    USED=OUT$USED
    DIST=OUT$DIST
    ALL_VEC=OUT$VEC
    USED_NAME=OUT$USED_NAME
    USED_CENTER_VEC=OUT$CENTER_VEC[USED,]
    USED_DIST=OUT$DIST[which(rownames(DIST) %in% USED_NAME),which(rownames(DIST) %in% USED_NAME)]
    OUT$SCORE=OUT$CENTER_VALUE[OUT$USED]
    SCORE=OUT$CENTER.VALUE[USED]
    ########################################
   
    #####################
    one=min(dist(USED_CENTER_VEC)) 
    #####################
    ###################
    .norm_one <-function(x,one=1){
        one=one
        if(var(x)!=0){
            x=x/sqrt(sum(x^2)) * one }
        return(x)
        }
    DIV=1/P
    ##########################
    N.SCORE=.norm1(SCORE)
    SCORE.COL=.vcol(N.SCORE, c(0,0.5,1),c('#009FFF','#FFF200','#ec2F4B'))
    
    A1_VEC=c()
    A2_VEC=c()
    #####################
    
    #####################
    i=1
    while(i<=length(USED)){
        this_p1_loc=USED_CENTER_VEC[i,]
        
        vector_list=cbind(USED_CENTER_VEC[,1]-this_p1_loc[1],USED_CENTER_VEC[,2]-this_p1_loc[2])
        vector_list_norm=t(apply(vector_list,1,.norm_one, one))
        
        vector_weight_1= DIV^-(rank(USED_DIST[i,])-1)   
        vector_weight_2= SCORE[i]-SCORE                  

        vector_weight_1=vector_weight_1 /sum(vector_weight_1)
        vector_weight_2=vector_weight_2 /sum(abs(vector_weight_2))
        
        vector_weight = vector_weight_1 * vector_weight_2        
        vector_weight = vector_weight/sum(abs(vector_weight))      
        final_vec=t(vector_list_norm) %*% vector_weight
        
        this_p2_loc=c(this_p1_loc[1]+final_vec[1],this_p1_loc[2]+final_vec[2])
            
        A1_VEC=cbind(A1_VEC,this_p1_loc)
        A2_VEC=cbind(A2_VEC,this_p2_loc)
        i=i+1
        }   
    ###############################
   
    #################################
    A1_VEC=t(A1_VEC)
    A2_VEC=t(A2_VEC)
    REV_VEC=A1_VEC-A2_VEC
    REV_VEC_CLOCK=.cart2clock(x=REV_VEC[,1],y=REV_VEC[,2])
    REV_VEC_CLOCK_UP=REV_VEC_CLOCK[,2]+30
    REV_VEC_CLOCK_DW=REV_VEC_CLOCK[,2]-30
    
    A2_VEC_UP=A2_VEC+.clock2cart(rho=REV_VEC_CLOCK[,1]*AL,phi=REV_VEC_CLOCK_UP)
    A2_VEC_DW=A2_VEC+.clock2cart(rho=REV_VEC_CLOCK[,1]*AL,phi=REV_VEC_CLOCK_DW)

    #################################
    if(SHOW==TRUE){
        if(BD==TRUE){
            plot(ALL_VEC,col=COL,pch=16,cex=CEX, 
                  xlim=c(min(ALL_VEC[,1])-one, max(ALL_VEC[,1])+one ),
                 ylim=c(min(ALL_VEC[,2])-one, max(ALL_VEC[,2])+one ) 
                )
            }else{
            plot(ALL_VEC,col=COL,pch=16,cex=CEX, 
                 xlim=c(min(ALL_VEC[,1])-one, max(ALL_VEC[,1])+one ),
                 ylim=c(min(ALL_VEC[,2])-one, max(ALL_VEC[,2])+one ) , yaxt="n", axes=F
                )
            }
        segments(x0=A1_VEC[,1],x1=A2_VEC[,1],y0=A1_VEC[,2],y1=A2_VEC[,2],lwd=LWD,col=AC)
        segments(x0=A2_VEC[,1],x1=A2_VEC_UP[,1],y0=A2_VEC[,2],y1=A2_VEC_UP[,2],lwd=LWD,col=AC)
        segments(x0=A2_VEC[,1],x1=A2_VEC_DW[,1],y0=A2_VEC[,2],y1=A2_VEC_DW[,2],lwd=LWD,col=AC) 

        }



    ######
    OUT$A1_VEC=A1_VEC
    OUT$A2_VEC=A2_VEC
    OUT$A2_VEC_UP=A2_VEC_UP
    OUT$A2_VEC_DW=A2_VEC_DW
    
    ###########
    return(OUT)
    }

 


#################################

fitdevo.field<-function(DP, VEC,COL=NULL, N=25, CUT=1, P=0.9, CEX=0.5, LWD=1.5, AL=0.4,SHOW=TRUE){
    DP=DP
    VEC=VEC
    TOSHOW=SHOW
    N=N
    CUT=CUT
    P=P
    COL=COL
    CEX=CEX
    LWD=LWD
    AL=AL
    #########################################
    if(length(DP)!=nrow(VEC)){
         print('DP and VEC are not matched !!!')
         return(0)
         }
    ################################################################
    this_field=.buildGrid(VEC, N=N, SHOW=FALSE)
    this_field=.buildNet(this_field, CUT=CUT, SHOW=FALSE)
    this_field$VALUE=DP
    this_field=.gridValue(this_field,SHOW=FALSE)
    if(is.null(COL)){COL=this_field$POINT.CENTER.COL}else{COL=COL}
    this_field=.drawArrow(this_field, P=P, SHOW=TOSHOW, COL=COL, CEX=CEX, LWD=LWD, AL=AL)
    return(this_field)
    }








#########################
# 2022.7.15
##########################


.simple_match <-function(exp_sc_mat1, exp_sc_mat2, FILL=FALSE){
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
    return(OUT)
    }


.norm1<-function(x){
    if(min(x)!=max(x)){
        y=(x-min(x))/(max(x)-min(x))
       }else{
        y=rep(0,length(x))
       }
    return(y)
    }



comdevo<-function(MAT, REF, DP=NULL,  PCNUM=5, NORM=TRUE, SEED=123, MAXDP=1){
    MAT=MAT
    DP=DP
    REF=REF
    NORM=NORM
    PCNUM=PCNUM
    SEED=SEED
    MAXDP=MAXDP
    ###########################################
    if(is.null(DP)){DP=rep(0,ncol(MAT))}
    ############################################
    library(Seurat)
    library(qlcMatrix)
    ###########################################
    REF.UP=.toUpper(REF)
    ###########################################
    #############################################
    # Create Seurat object
    pbmc=CreateSeuratObject(counts = MAT, min.cells = 0, min.features = 0, project = "ALL")
    pbmc$dp=DP

    #############################################
    # Normalization
    if(NORM==TRUE){  
        pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000) 
        }

    #################################################
    NMAT=pbmc@assays$RNA@data
    NCOL=ncol(NMAT)
    NMAT.UP=.toUpper(NMAT)
    ####################################################
     
    MATCH=.simple_match(REF.UP, NMAT.UP)
    D1=MATCH$exp_sc_mat1
    D2=MATCH$exp_sc_mat2

    PCC.MAT=Matrix(corSparse(Matrix(D1),Matrix(D2)))
    rownames(PCC.MAT)=colnames(D1)
    colnames(PCC.MAT)=colnames(D2)
    INPUT=(PCC.MAT+1) *50

    ####################################################
    pbmc[["PCC"]] <- CreateAssayObject(counts = INPUT )
    DefaultAssay(pbmc)='PCC'
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000) 
    all.genes=rownames(pbmc)
    pbmc <- ScaleData(object = pbmc, features = all.genes)
    pbmc <- RunPCA(object = pbmc, seed.use=SEED, npcs=PCNUM, features = all.genes, ndims.print=1,nfeatures.print=1)

    #################################
    this_pca=pbmc@reductions$pca@cell.embeddings

    this_weight=MAXDP-pbmc$dp
    this_pca=this_pca * this_weight

    #############################################

    ##############################################
    rownames(this_pca)=rownames(pbmc@reductions$pca@cell.embeddings)
    colnames(this_pca)=colnames(pbmc@reductions$pca@cell.embeddings)

    #############################################
    pbmc@reductions$pca@cell.embeddings=this_pca

    ##################################
    pbmc <- RunUMAP(pbmc, dims = 1:PCNUM, seed.use = SEED, n.components=2)  
    
    #DimPlot(pbmc,label=TRUE)+NoLegend()
    #FeaturePlot(pbmc,features=c('dp'))

    ##################################
    #DefaultAssay(pbmc)='RNA'
    return(pbmc)

    }


