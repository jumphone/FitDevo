<img src="https://github.com/jumphone/FitDevo/blob/main/img/logo.png?raw=true" width="250">

**FitDevo: accurate inference of single-cell developmental potential using sample-specific gene weight, Briefings in Bioinformatics, 2022, in press** 

This tool is designed for inferring the developmental potential (DP) of cells in scRNA-seq data

# Updates:

**2022.07.01, v1.1.0 </br></br> New features!** Users can use "fitdevo.field" to build developmental potential field (DPF) and draw arrows.

**2022.06.30, v1.0.1 </br></br> Paper version.** The details of this version is described in our BIB (2022) paper.

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
    BGW: binarized gene weight (BGW) ('https://github.com/jumphone/FitDevo/blob/main/BGW.rds?raw=true')
    NORM: whether to run "LogNormalize" in Seurat
    PCNUM: number of PCs used to calculate sample-specific gene weight (SSGW)

## Output:

    A vector of inferred DP
    
--------------------------------------------------------------------------------------------------------------------

## Demo 1 | Infer developmental potential (DP) using expression matrix of scRNA-seq data

    # R 4.0.3 
    
    # Step 1. Load FitDevo 
    source('https://github.com/jumphone/FitDevo/blob/main/fitdevo.R?raw=true')
    
    # Step 2. Load data (the 1st sample in the testing dataset)
    data1 = readRDS('1.rds')
    MAT=data1$mat
    CorrectDP=data1$tag
    
    # Step 3. Load BGW
    BGW=readRDS(url('https://github.com/jumphone/FitDevo/blob/main/BGW.rds?raw=true'))
    
    # Step 4. Run FitDevo
    DP=fitdevo(MAT=MAT, BGW=BGW, NORM=TRUE, PCNUM=50)
    
    # Step 5. Evaluate the performance of FitDevo
    cor(DP, CorrectDP, method='spearman')  # 0.7980606

--------------------------------------------------------------------------------------------------------------------

## Demo 2 | Build developmental potential field (DPF) and draw arrows ( fitdevo >= 1.1.0 )
    
Users should provide embedding coordinates (e.g. tSNE, UMAP, PAGA, etc.).

This demo is based on a "seurat.object" with normalized expression matrix and UMAP.

To generate seurat.object, please refer to: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html
    
    source('https://github.com/jumphone/FitDevo/blob/main/fitdevo.R?raw=true')
    
    # Step 1. Prepare input files.
    MAT=as.matrix(seurat.object[['RNA']]@data)
    VEC=seurat.object@reductions$umap@cell.embeddings
    BGW=readRDS(url('https://github.com/jumphone/FitDevo/blob/main/BGW.rds?raw=true'))
    
    # Step 2. Infer developmental potential
    DP=fitdevo(MAT, BGW, NORM=FALSE, PCNUM=50)

    # Step 3. Build developmental potential field (DPF) and draw arrows
    FIELD=fitdevo.field(DP=DP, VEC=VEC, SHOW=TRUE)


<img src="https://github.com/jumphone/FitDevo/blob/main/img/f01_demo2_fitdevo.field.png?raw=true" width="300">

--------------------------------------------------------------------------------------------------------------------
    


