<img src="https://gitee.com/jumphone/public/raw/master/fitdevo_logo.png" width="250">

FitDevo: accurate inference of single-cell developmental potential using sample-specific gene weight

This tool is designed for inferring developmental potential (DP) using scRNA-seq data


# Training & Testing Datasets:

Training dataset (n=17), https://sourceforge.net/projects/fitdevo/files/training/

Testing dataset (n=28), https://sourceforge.net/projects/fitdevo/files/testing/

Each sample is saved in a "RDS" file. 

Users can use R to load the "RDS" file.
    
    # R code
    data1 = readRDS('1.rds')
    # data1$mat is the expression matrix
    # data1$tag is the reverse order of timepoint label (higher value indicates higher developmental potential).
    

# Requirements:

    R: 4.0.0+
    Seurat: 4.0+
    
    
# Usage:

    fitdevo( MAT=MAT, BGW=BGW, NORM=TRUE, PCNUM=50 )

## Input:

    MAT: expression matrix
    BGW: binarized gene weight (BGW) ('https://gitee.com/jumphone/public/raw/master/BGW.rds' or 'https://github.com/jumphone/FitDevo/blob/main/BGW.rds?raw=true')
    NORM: whether to run "LogNormalize" in Seurat
    PCNUM: number of PCs used to calculate sample-specific gene weight (SSGW)

## Output:

    A column vector of inferred DP

## Demo:

    # R 4.0.3 
    
    # Step 1. Load FitDevo 
    source('https://gitee.com/jumphone/public/raw/master/fitdevo.R') 
    # or source('https://github.com/jumphone/FitDevo/blob/main/fitdevo.R?raw=true')
    
    # Step 2. Load data (the 1st sample in testing dataset)
    data1 = readRDS('1.rds')
    MAT=data1$mat
    CorrectDP=data1$tag
    
    # Step 3. Load BGW
    BGW=readRDS(url('https://gitee.com/jumphone/public/raw/master/BGW.rds')) 
    # or BGW=readRDS(url('https://github.com/jumphone/FitDevo/blob/main/BGW.rds?raw=true'))
    
    # Step 4. Run FitDevo
    DP=fitdevo(MAT=MAT, BGW=BGW, NORM=TRUE, PCNUM=50)
    
    # Step 5. Evaluate the performance of FitDevo
    cor(DP, CorrectDP, method='spearman')  # 0.7980606

