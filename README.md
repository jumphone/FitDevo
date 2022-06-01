<img src="https://gitee.com/jumphone/public/raw/master/fitdevo_logo.png" width="250">

FitDevo: accurate inference of single-cell developmental potential using sample-specific gene weight

This tool is designed for inferring the developmental potential (DP) of cells in scRNA-seq data


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
    
The input of FitDevo includes two files: BGW file and expression matrix. The BGW file is provided by us, while the expression matrix is provided by users. The row and column names of the expression matrix are gene and cell names, respectively. FitDevo can help users to normalize the raw read count by setting the “NORM” parameter to “TRUE”, or users can use “LogNormalize” function in Seurat to conduct normalization. The output of FitDevo is a vector containing the inferred DP of all cells.

## Input:

    MAT: expression matrix
    BGW: binarized gene weight (BGW) ('https://gitee.com/jumphone/public/raw/master/BGW.rds' or 'https://github.com/jumphone/FitDevo/blob/main/BGW.rds?raw=true')
    NORM: whether to run "LogNormalize" in Seurat
    PCNUM: number of PCs used to calculate sample-specific gene weight (SSGW)

## Output:

    A column vector of inferred DP

## Demo 1: use FitDevo alone:

    # R 4.0.3 
    
    # Step 1. Load FitDevo 
    source('https://gitee.com/jumphone/public/raw/master/fitdevo.R') 
    # or source('https://github.com/jumphone/FitDevo/blob/main/fitdevo.R?raw=true')
    
    # Step 2. Load data (the 1st sample in the testing dataset)
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




## Demo 2: combine FitDevo with trajectory-based method (Vector):

Vector: https://github.com/jumphone/Vector

    # R 4.0.3 
    
    # Step 1. Load FitDevo 
    source('https://gitee.com/jumphone/public/raw/master/fitdevo.R') 
    # or source('https://github.com/jumphone/FitDevo/blob/main/fitdevo.R?raw=true')
    
    # Step 2. Load data
    seuratData = readRDS('brain.rds')
    
    # Step 3. Load BGW
    BGW=readRDS(url('https://gitee.com/jumphone/public/raw/master/BGW.rds')) 
    # or BGW=readRDS(url('https://github.com/jumphone/FitDevo/blob/main/BGW.rds?raw=true'))
    
    # Step 4. Run FitDevo
    MAT=seuratData[['RNA']]@counts
    DP=fitdevo(MAT=MAT, BGW=BGW, NORM=TRUE, PCNUM=50)
    
    # Step 5. Load Vector
    source('https://gitee.com/jumphone/Vector/raw/master/Vector.R')
    # or source('https://raw.githubusercontent.com/jumphone/Vector/master/Vector.R') 
    
    # Step 6. Combine FitDevo & Vector
    VEC = seuratData@reductions$umap@cell.embeddings
    rownames(VEC) = colnames(seuratData)
    PCA = seuratData@reductions$pca@cell.embeddings
    
    OUT=vector.buildGrid(VEC, N=30,SHOW=FALSE)
    OUT=vector.buildNet(OUT, CUT=1, SHOW=FALSE)
    OUT$VALUE=DP
    OUT=vector.gridValue(OUT,SHOW=FALSE)
    OUT=vector.autoCenter(OUT,UP=0.9,SHOW=FALSE)
    OUT=vector.drawArrow(OUT,P=0.9,SHOW=TRUE, COL=OUT$COL, SHOW.SUMMIT=TRUE)

    
<img src="https://gitee.com/jumphone/public/raw/master/fitdevo_vector_tra.jpg" width="250">    
<img src="https://gitee.com/jumphone/public/raw/master/fitdevo_type.jpg" width="250">    
    
    
