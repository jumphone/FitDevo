
# FitDevo: accurate inference of single-cell developmental potential using sample-specific gene weight

FitDevo is designed for inferring developmental potential (DP) using scRNA-seq data










# Datasets

Training dataset (n=17):

https://sourceforge.net/projects/fitdevo/files/training/

Testing dataset (n=28):

https://sourceforge.net/projects/fitdevo/files/testing/


Each sample is saved in a "RDS" file. Users can use R to load the "RDS" file.
    
    # R code
    data1 = readRDS('1.rds')
    # data1$mat is the expression matrix
    # data1$tag is the differentiation label. Higher value indicates higher developmental potential.
    
    
    
    
