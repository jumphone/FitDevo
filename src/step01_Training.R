

# /home/toolkit/tools/R4.0.3/bin/R


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


