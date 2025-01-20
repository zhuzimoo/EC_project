import scanpy as sc
import anndata as ad
import os
import scanpy.external as sce

index_list = [i for i in range(0, 17)]
ori = sc.read_h5ad("./raw_used_tissues/disco_adipose_v01.h5ad")
adata = sc.AnnData(X=ori.X, obs=ori.obs, var=ori.var)
adata = adata[(adata.obs.ct == 'Arterial EC') |(adata.obs.ct == 'Venous EC') | (adata.obs.ct == 'Capillary EC') | (adata.obs.ct == 'Lymphatic EC')]
adata = adata[adata.obs.disease == "NA"]
adata.obs["tissue_comb"] = "adipose"
adata.obs["index_num"] = str(index_list[0])
if "project_id" in adata.obs_keys():
    adata.obs["project_comb"] = adata.obs["project_id"]
else:
    adata.obs["project_comb"] = adata.obs["projectId"]
print(adata.obs_keys())

results_folder = "./raw_used_tissues"
j = 1
directory = os.fsencode(results_folder)
for filename in sorted(os.listdir(directory)):
    filename = os.fsdecode(filename)
    if "adipose" in filename:
        continue
    path = results_folder + "/" + filename
    tissue = sc.read_h5ad(path)
    tissue_name = filename.split('_')[1]
    adata_new = sc.AnnData(X=tissue.X, obs=tissue.obs, var=tissue.var)
    adata_new.obs["tissue_comb"] = tissue_name
    adata_new.obs["index_num"] = str(index_list[j])
    if "disease" in list(adata_new.obs_names):
        adata_new = adata_new[adata_new.obs.disease == "NA"]
    adata_new = adata_new[(adata_new.obs.ct == 'Arterial EC') |(adata_new.obs.ct == 'Venous EC') | (adata_new.obs.ct == 'Capillary EC') | (adata_new.obs.ct == 'Lymphatic EC')]
    print(adata_new.shape)
    if "project_id" in adata.obs_keys():
        adata.obs["project_comb"] = adata.obs["project_id"]
    else:
        adata.obs["project_comb"] = adata.obs["projectId"]
    adata = ad.concat([adata, adata_new], join="outer")
    print(adata.shape)
    j += 1
sc.set_figure_params(dpi=600)
adata.write('write/comb_ori_4ec.h5ad')

# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
adata.write('write/comb_norm_4ec.h5ad')

# Batch removal
sce.pp.harmony_integrate(adata, 'project_comb')
sc.pp.neighbors(adata, use_rep="X_pca_harmony")
sc.tl.umap(adata)
adata.write('write/comb_br_4ec.h5ad')


