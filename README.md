<img src="https://fzhang.bioinfo-lab.com/img/tools/logo_fitdevo.png" width="250">

**FitDevo: accurate inference of single-cell developmental potential using sample-specific gene weight**, ***Briefings in Bioinformatics*, 2022**

**Paper Link:** [https://doi.org/10.1093/bib/bbac293](https://doi.org/10.1093/bib/bbac293) | [Supplementary Files](https://github.com/jumphone/FitDevo/tree/main/sup)

This tool is designed for inferring the developmental potential (DP) of cells in scRNA-seq data


</br>

# Updates:

**2022.07.13, v1.2.0 - Higher speed!** Improve the computational efficiency by using "qlcMatrix".

**2022.07.01, v1.1.0 - New features!**  Users can use "fitdevo.field" to build developmental potential field (DPF) and draw arrows.

**2022.06.30, v1.0.1 - Paper version.** The details of this version is described in our BIB (2022) paper.


</br>


# Content:

* [Training & Testing Datasets](#training--testing-datasets)
* [Requirements](#requirements)
* [Usage](#usage)

### Demos:

* [Demo 1 - Infer developmental potential (DP) using expression matrix of scRNA-seq data](#demo-1---infer-developmental-potential-dp-using-expression-matrix-of-scrna-seq-data)
* [Demo 2 - Build developmental potential field (DPF) and draw arrows](#demo-2---build-developmental-potential-field-dpf-and-draw-arrows--fitdevo--110-)

</br>
</br>
</br>



# Training & Testing Datasets:

Training dataset (n=17), https://sourceforge.net/projects/fitdevo/files/training/

Testing dataset (n=28), https://sourceforge.net/projects/fitdevo/files/testing/

Each sample is saved in a "RDS" file. 

Users can use R to load the "RDS" file.
    
    # R code
    data1 = readRDS('1.rds')
    # data1$mat is the expression matrix
    # data1$tag is the reverse order of timepoint label (higher value indicates higher developmental potential).
    
</br>

# Requirements:

    R: 4.0.0+
    Seurat: v4 (the original fitdevo only supports Seurat v4)
    qlcMatrix: 0.9.7

FitDevo for Seurat v5: https://github.com/jumphone/FitDevo/blob/main/fitdevo_seurat5.R

R: https://www.r-project.org/

Seurat: https://satijalab.org/seurat/articles/install.html

qlcMatrix: https://cran.r-project.org/web/packages/qlcMatrix/index.html

    install.packages("qlcMatrix")

Seurat v4:

    remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
    remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

Seurat v5:

    install.packages("Seurat") 
    


</br>


# Usage:
    
    source('https://github.com/jumphone/FitDevo/blob/main/fitdevo.R?raw=true') 

    fitdevo( MAT=MAT, BGW=BGW, NORM=TRUE, PCNUM=50 )
    
The input of FitDevo includes two files: a BGW list and an expression matrix. The BGW list is provided by us, while the expression matrix is provided by users (should not be scaled). The row and column names of the expression matrix are genes and cell names, respectively. FitDevo can help users to normalize the raw read count by setting the “NORM” parameter to “TRUE”, or users can use “LogNormalize” function in Seurat to conduct normalization. The output of FitDevo is a vector containing the inferred DP of all cells.

### Input:

    MAT: expression matrix
    BGW: binarized gene weight (BGW)
    
    BGW=readRDS(url('https://github.com/jumphone/FitDevo/blob/main/BGW.rds?raw=true'))
    
    NORM: whether to conduct normalization
    PCNUM: number of PCs used to calculate sample-specific gene weight (SSGW)

### Output:

    A vector of inferred DP
    
</br>

[Click back to the top](#)

</br>
</br>



# Demos:

</br>

## Demo 1 - Infer developmental potential (DP) using expression matrix of scRNA-seq data

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


</br>

[Click back to the top](#)

</br>
</br>


## Demo 2 - Build developmental potential field (DPF) and draw arrows ( latest version of FitDevo )

### Please install "igraph" and "stringr" before using "fitdevo.field"

    install.packages('igraph')
    install.packages('stringr')

Users should provide the embedding coordinates (e.g. tSNE, UMAP, PAGA, etc.). This demo is based on a "seurat.object" with a normalized expression matrix and an UMAP. To generate seurat.object, please refer to: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

The function named "fitdevo.field" partially follows the idea of another tool named VECTOR. If you are using this function in your work, please also cite: Unsupervised Inference of Developmental Directions for Single Cells Using VECTOR, Cell Reports, 2020. {[code](https://github.com/jumphone/Vector), [paper](https://doi.org/10.1016/j.celrep.2020.108069)}

</br>

    source('https://github.com/jumphone/FitDevo/blob/main/fitdevo.R?raw=true')
    
    # Step 1. Prepare input files.
    MAT=as.matrix(seurat.object[['RNA']]@data)
    VEC=seurat.object@reductions$umap@cell.embeddings
    BGW=readRDS(url('https://github.com/jumphone/FitDevo/blob/main/BGW.rds?raw=true'))
    
    # Step 2. Infer developmental potential
    DP=fitdevo(MAT, BGW, NORM=FALSE, PCNUM=50)

    # Step 3. Build developmental potential field (DPF) and draw arrows
    FIELD=fitdevo.field(DP=DP, VEC=VEC, SHOW=TRUE)
    
    # Coler scale
    plot(x=FIELD$VALUE,y=rep(1,length(FIELD$VALUE)),col=FIELD$POINT.COL,type='h',lwd=2,ylim=c(0,1))
    
<p float="left">
<img src="https://raw.githubusercontent.com/jumphone/FitDevo/main/img/f01_demo2_fitdevo.field.png" width="250">
</p>

This figure is generated by using the scRNA-seq data of mouse dentate gyrus ([PMID: 29335606](https://www.nature.com/articles/s41593-017-0056-2)).

</br>

[Click back to the top](#)

</br>
</br>



<img src="https://fzhang.bioinfo-lab.com/img/white.png" height="50">

-------------------------------------------------------------------------------------------------------------------

<img src="https://fzhang.bioinfo-lab.com/img/panda_happy_logo.png" height='150'>

#### More tools & studies: https://fzhang.bioinfo-lab.com/
