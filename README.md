
FitDevo: accurate inference of single-cell developmental potential using sample-specific gene weight

FitDevo is designed for inferring developmental potential (DP) using scRNA-seq data


# Training & Testing Datasets:

Training dataset (n=17):

https://sourceforge.net/projects/fitdevo/files/training/

Testing dataset (n=28):

https://sourceforge.net/projects/fitdevo/files/testing/


Each sample is saved in a "RDS" file. Users can use R to load the "RDS" file.
    
    # R code
    data1 = readRDS('1.rds')
    # data1$mat is the expression matrix
    # data1$tag is the differentiation label. Higher value indicates higher developmental potential.
    

# Requirements:

    R: 4.0.0+
    Seurat: 4.0+
    
    
# Installation:
    
    # Directly load the code in R
    source('https://raw.githubusercontent.com/jumphone/fitdevo/main/fitdevo.R') # or source('')
    
# Usage:

    #/home/toolkit/tools/R4.0.3/bin/R
    source('https://gitee.com/jumphone/public/raw/master/fitdevo.R')
    
    # Load the first sample in testing dataset
    data1 = readRDS('1.rds')
    MAT=data1$mat
    CorrectDP=data1$tag
    
    # Load BGW
    BGW=readRDS(url('https://gitee.com/jumphone/public/raw/master/BGW.rds'))
    
    # Run FitDevo
    DP=fitdevo(MAT=MAT, BGW=BGW, NORM=TRUE, PCNUM=50)
    
    # Evaluate the performance of FitDevo
    cor(DP,CorrectDP, method='spearman')
    # 0.7980606

