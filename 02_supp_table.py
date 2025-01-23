import scanpy as sc
import anndata as ad
import os
import pandas as pd


comb_table = pd.DataFrame()

# load in raw scRNA-seq data
results_folder = "./raw_used_tissues"
directory = os.fsencode(results_folder)

for filename in sorted(os.listdir(directory)):
    filename = os.fsdecode(filename)
    path = results_folder + "/" + filename
    tissue_name = filename.split('_')[1]
    ori = sc.read_h5ad(path)
    adata = sc.AnnData(X=ori.X, obs=ori.obs, var=ori.var)
    adata.obs["tissue"] = tissue_name
    if "project_id" in adata.obs.columns.tolist():
        table = adata.obs[["tissue", "project_id"]]
        table["cell_num"] = table.groupby("project_id").transform("count")
        table = table.drop_duplicates(subset=["project_id", "tissue"])
        table = table.reset_index(drop=True)
    else:
        table = adata.obs[["tissue", "projectId"]]
        table["cell_num"] = table.groupby("projectId").transform("count")
        table = table.drop_duplicates(subset=["projectId", "tissue"])
        table = table.reset_index(drop=True)
        table = table.rename(columns={"projectId": "project_id"})
    print(table)
    comb_table = pd.concat([table, comb_table])
print(comb_table)
comb_table.to_csv("project_id_info_new.csv", index=False)
