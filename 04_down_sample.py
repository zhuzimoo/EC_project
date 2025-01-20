import scanpy as sc
import anndata as ad
import os

tissue_name = "adipose"
export_file = "./raw_10000/" + tissue_name + "_10000.h5ad"
ori = sc.read_h5ad("./DISCO/raw_used_tissues/disco_thymus_v01.h5ad")
adata_X = ori.X
adata_obs = ori.obs
adata_var = ori.var
del ori
adata = sc.AnnData(X=adata_X, obs=adata_obs, var=adata_var)
print(adata)
if "sample_type" in adata.obs.columns.tolist():
   adata = adata[adata.obs.disease == "NA"]
   adata = adata[(adata.obs.sample_type == "Normal") | (adata.obs.sample_type == "normal")| (adata.obs.sample_type == "Adjacent normal") | (adata.obs.sample_type == "adjacent normal")]
else:
   adata = adata[adata.obs.disease == "NA"]
   adata = adata[(adata.obs.sampleType == "Normal") | (adata.obs.sampleType == "normal")| (adata.obs.sampleType == "Adjacent normal") | (adata.obs.sampleType == "adjacent normal")]

adata.obs["tissue_comb"] = tissue_name
sc.pp.subsample(adata, n_obs=10000, random_state=123)
adata.write_h5ad(export_file)

