# scripts for scRNA-seq data analysis
### Dependency and version control
- Python
    - processing raw scRNA-seq data: Python==3.8
        - packages: Scanpy (version 1.8.2), Anndata (version 0.9.2)
    - downstream analysis: Python==3.10
        - packages: Scipy (version 1.13.0), scikit-learn (version 1.4.2), statsmodels (version 0.14.2)
- R==4.4
    - packages: topicmodels (version 0.2-16), tidytext (version 0.4.2), tidy (version 1.3.1), dplyr (version 1.1.4), pheatmap (version 1.0.12), ggplot2 (version 3.5.1)

### Installation
The packages and softwares can be installed by the following approaches:
- Anaconda. For example, scanpy can by installed from anaconda at https://anaconda.org/bioconda/scanpy
- Python packages can be installed from `pip install` command lines
- R packages can be installed from `install.packages('<package_name>')` command lines
- [MEBOCOST](https://github.com/kaifuchenlab/MEBOCOST), [LIANA](https://liana-py.readthedocs.io/en/latest/index.html)

### Downloading the data
The raw scRNA-seq data could be found at [DISCO](https://www.immunesinglecell.org/). The data and corresponding version used is listed as follows: 

adipose (version 1.0), bladder (version 1.1), breast (version 2.1), gut (version 1.0), heart (version 1.0), intestine (version 1.0), kidney (version 1.0), liver (version 2.0), lung (version 2.0), ovary (version 1.0), skeletal muscle (version 1.0), skin (version 1.0), stomach (version 1.0), testis (version 1.0), thymus (version 1.0) [all in H5AD format]

### Table of content
- 01_combine_adata.py

    The python script for combining anndata from all dataset downloaded, normalization, and batch effect removal.

- 02_supp_table.py

    The python script for generating the supplementary table for project and cell information stored in the scRNA-seq data.

- 03_pca_analysis.py

    The python script for PCA analysis and hierarchical clustering using the 50 harmony-corrected PCs.

- 04_down_sample.py

    The python script for raw scRNA-seq data downsampling for each individual data set.

- 05_MEBOCOST_run.py

    The python script for running MEBOCOST to infer metabolite-mediated cell-cell communication (CCC) events in each data set. 

- 06_LIANA_run.py

    The python script for running LIANA to infer protein-mediated CCC events in each data set.

- 07_commu_analysis.py

    The python script for analysis and comparison in the inferred communication events under two contexts. 

- 08_LDA_matrix_MEBOCOST.py

    The python script for generating input matrix for LDA analysis using the data inferred by MEBOCOST. 

- 09_LDA_matrix_LIANA.py

    The python script for generating input matrix for LDA analysis using the data inferred by LIANA. 

- 10_LDA.R

    The R script for LDA analysis using topic modeling. 

- 11_LDA_analysis_MEBOCOST.py

    The python script for downstream analysis based on LDA output matrix with metabolite-mediated CCC events.

- 12_LDA_analysis_LIANA.py

    The python script for downstream analysis based on LDA output matrix with protein-mediated CCC events. 



