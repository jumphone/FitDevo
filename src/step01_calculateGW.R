

# /home/toolkit/tools/R4.0.3/bin/R

library(Seurat)
.normData<-function(mat){
    mat=mat
    mat=mat[which(rownames(mat)!=''),]
    pbmc=CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0, project = "ALL")
    pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    this_out=as.matrix(pbmc@assays$RNA@data)
    return(this_out)
    }


.simple_combine_NA <- function(exp_sc_mat1, exp_sc_mat2, FILL=FALSE){
    FILL=FILL
    exp_sc_mat=exp_sc_mat1
    exp_ref_mat=exp_sc_mat2
    ##############################################
    if(FILL==TRUE){
        gene1=rownames(exp_sc_mat)
        gene2=rownames(exp_ref_mat)
        gene12=gene2[which(!gene2 %in% gene1)]
        gene21=gene1[which(!gene1 %in% gene2)]
        exp_sc_mat_add=matrix(NA,ncol=ncol(exp_sc_mat),nrow=length(gene12))
        rownames(exp_sc_mat_add)=gene12
        colnames(exp_sc_mat_add)=colnames(exp_sc_mat)
        exp_ref_mat_add=matrix(NA,ncol=ncol(exp_ref_mat),nrow=length(gene21))
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



##########################################################################################
# Get Protein Coding Genes
##########################################################################################

SYM_BED=read.table('/home/zhangfeng/project/DEV/data/annotation/Homo_sapiens.GRCh37.75.chr.pc.gene.SYM.bed',header=F,sep='\t')
PCGENE=SYM_BED[,4]

head(PCGENE)
# "OR4F5" "AL627309.1" "OR4F29" "OR4F16" "AL669831.1" "AL645608.2"



##########################################################################################
# Get "gene expression" & "reverse time order" correlation matrix of 17 samples collected by the study of CytoTRACE
##########################################################################################

setwd('/home/zhangfeng/project/DEV/data/CytoTRACE')

TMP=rep(0,length(PCGENE))
names(TMP)=PCGENE
TMP=cbind(TMP,TMP)
TMP=.toUpper(TMP)


GCC=TMP

i=1
while(i<=17){
        ID=i
        this_file=paste0(ID,'.rds')
        this_data=readRDS(this_file)
        this_mat=this_data$mat
        ##############
        this_mat= .normData(this_mat)
        ################################################
        #######################################
        this_mat=.toUpper(this_mat)
        #############
        this_tag= this_data$tag
        ####################################
        used_index=which(!is.na(this_tag))
        this_mat=this_mat[,used_index]
        this_tag=this_tag[used_index]
        ############################################################################
        #####################################
        this_mat=as.matrix(this_mat)
        #####################################
        this_mat=this_mat[which(rownames(this_mat)!=''),]
        #####################################
        this_exp=this_mat[which(rownames(this_mat) %in% rownames(TMP)),]
        #############################
        this_c=cor(t(this_exp), this_tag)
        ############################
        NGCC=ncol(GCC)
        GCC=.simple_combine_NA(GCC,cbind(this_c,this_c), FILL=TRUE)$combine[,c(1:(NGCC+1))]     
        print(i)
        i=i+1
       }

GCC=GCC[,c(3:ncol(GCC))]

MMM=apply(GCC,1,mean,na.rm=TRUE)
################
TMP=GCC
TMP[which(!is.na(GCC))]=0
TMP[which(is.na(GCC))]=1
NAN=apply(TMP,1,sum)
##############
GCC=GCC[which(NAN<=8 & MMM!=0),]
####################

S.GCC=t(apply(GCC,1,scale.rmNA, center=FALSE)) # dividing Root-Mean-Square (RMS)
rownames(S.GCC)=rownames(GCC)



#######################################################################

GW=apply(S.GCC,1,mean,na.rm=TRUE) # Gene Weight (GW)

#######################################################################

# Binarized Gene Weight (BGW)
BGW=GW
BGW[which(GW>0)]=1
BGW[which(GW<=0)]=0

########################################################################


