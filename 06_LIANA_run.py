import liana as li
import scanpy as sc
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean
import matplotlib.pyplot as plt

tissue_name = "thymus"
adata = sc.read_h5ad(f"./raw_10000/{tissue_name}_10000.h5ad")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# run cellphoned
cellphonedb(adata, groupby='ct', expr_prop=0.1, resource_name='consensus', verbose=True, key_added='cpdb_res')
print(adata.uns['cpdb_res'].head())

adata.uns['cpdb_res'].to_csv(f"./comb_result/liana/{tissue_name}_cpdb_res.csv", index=False)

# run LIANA aggregate
li.mt.rank_aggregate(adata, groupby='ct', expr_prop=0.1, verbose=True)
adata.uns['liana_res'].to_csv(f"./comb_result/liana/{tissue_name}_liana_res.csv", index=False)

