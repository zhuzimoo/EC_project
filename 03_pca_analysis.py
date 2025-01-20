import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import pairwise_distances
import numpy as np

df = pd.read_csv("harmony_pc_tissue_ct.csv")
df.groupby(['Cell_type', 'Tissue_type']).mean()
pd.DataFrame(df.groupby(['Cell_type', 'Tissue_type']).mean()).to_csv("harmony_avg_pc.csv")

avg_pc = pd.read_csv("harmony_avg_pc.csv").reindex().sort_values("Tissue_type")

# Fig.1E
plt.figure(figsize=(8,6))
sns.set_style(style='white')
sns.scatterplot(x="PC1", y="PC2", data=avg_pc,
                alpha=0.8, hue="Tissue_type",  s=300, style="Cell_type",
                palette=["#1f77b4", "#aec7e8", "#ff7f0e", "#2ca02c", "#98df8a", "#ff9896", "#9467bd", "#8c564b",
                         "#c49c94", "#e377c2", "#7f7f7f", "#c7c7c7", "#bcbd22", "#17becf", "#9edae5"])

plt.xlabel('PC1')
plt.ylabel('PC2')
plt.legend(loc='upper left', bbox_to_anchor=(1.02, 1))
plt.title("PCA Analysis")
plt.tight_layout()
# plt.savefig()


new = avg_pc.copy()
new["tissue_ct"] = new["Tissue_type"] + "_" + new["Cell_type"]
new["tissue_ct"] = new["tissue_ct"].str.replace(" EC", "")

# Fig.1F
plt.figure(figsize=(5,3))
index = new[["tissue_ct"]]
index = index.to_numpy()
linkage_data = linkage(new.iloc[:, [2, 51]], method='ward', metric='euclidean')
dendrogram_data = dendrogram(linkage_data, labels=index)
plt.tight_layout()
# plt.savefig()

# Divide the major two clusters
cluster_info = {key: dendrogram_data[key] for key in ["ivl", "leaves_color_list"]}
cluster_info["ivl"]
cluster_key_list = []
for each in cluster_info["ivl"]:
    each = "".join(np.array_str(each)).replace("[", "").replace("]", "").replace("'", "")
    cluster_key_list.append(each)
cluster_num_list = cluster_info["leaves_color_list"]

cluster_df = pd.DataFrame({"cluster_key": cluster_key_list,
                           "cluster_num": cluster_num_list})
c1_list = cluster_df[cluster_df["cluster_num"] == "C1"]["cluster_key"].tolist()
c2_list = cluster_df[cluster_df["cluster_num"] == "C2"]["cluster_key"].tolist()

# Cluster 1
new_tissue_c1 = new.copy()
new_tissue_c1 = new_tissue_c1[new_tissue_c1["tissue_ct"].isin(c1_list)]
new_tissue_c1 = new_tissue_c1.sort_values("Cell_type")
tissue_distance_matrix_c1 = []
for each_tissue in new_tissue_c1.Tissue_type.unique().tolist():
    print(each_tissue)
    tissue_matix = new_tissue_c1[new_tissue_c1["Tissue_type"] == each_tissue]
    tissue_matix = tissue_matix.drop(columns=["Tissue_type", "tissue_ct"])
    d_matrix = pairwise_distances(tissue_matix.iloc[:, 1:], metric='euclidean')
    print(d_matrix)
    d_list = list(set(i for j in d_matrix for i in j))
    d_list.remove(0)
    print(d_list)
    tissue_distance_matrix_c1.append(d_list)
    print()
print([x for xs in tissue_distance_matrix_c1 for x in xs])

new_ct_c1 = new.copy()
new_ct_c1 = new_ct_c1[new_ct_c1["tissue_ct"].isin(c1_list)]
new_ct_c1 = new_ct_c1.sort_values("Tissue_type")
ct_distance_matrix_c1 = []
for each_ct in new_ct_c1.Cell_type.unique().tolist():
    print(each_ct)
    ct_matix = new_ct_c1[new_ct_c1["Cell_type"] == each_ct]
    ct_matix = ct_matix.drop(columns=["Cell_type", "tissue_ct"])
    print(ct_matix)
    d_matrix = pairwise_distances(ct_matix.iloc[:, 1:], metric='euclidean')
    print(d_matrix)
    d_list = list(set(i for j in d_matrix for i in j))
    d_list.remove(0)
    print(d_list)
    ct_distance_matrix_c1.append(d_list)
    print()
print([x for xs in ct_distance_matrix_c1 for x in xs])

# Fig.1G
plt.figure().set_figwidth(2.2)
sns.set(style="ticks")
boxplot_data = [[x for xs in ct_distance_matrix_c1 for x in xs], [x for xs in tissue_distance_matrix_c1 for x in xs]]
_ = plt.boxplot(boxplot_data, patch_artist = True, boxprops = dict(facecolor = "orange"), medianprops = dict(color = "red"),
                widths=(0.4, 0.4))
_ = plt.xticks([1, 2], ['Tissue type', 'EC type'])
_ = plt.ylabel("Distance")
_ = plt.title("cluster 1")
plt.tight_layout()
plt.ylim(0, 25)

# Cluster 2
new_tissue_c2 = new.copy()
new_tissue_c2 = new_tissue_c2[new_tissue_c2["tissue_ct"].isin(c2_list)]
new_tissue_c2 = new_tissue_c2.sort_values("Cell_type")
tissue_distance_matrix_c2 = []
for each_tissue in new_tissue_c2.Tissue_type.unique().tolist():
    print(each_tissue)
    tissue_matix = new_tissue_c2[new_tissue_c2["Tissue_type"] == each_tissue]
    tissue_matix = tissue_matix.drop(columns=["Tissue_type", "tissue_ct"])
    d_matrix = pairwise_distances(tissue_matix.iloc[:, 1:], metric='euclidean')
    print(d_matrix)
    d_list = list(set(i for j in d_matrix for i in j))
    d_list.remove(0)
    print(d_list)
    tissue_distance_matrix_c2.append(d_list)
    print()
print([x for xs in tissue_distance_matrix_c2 for x in xs])

new_ct_c2 = new.copy()
new_ct_c2 = new_ct_c2[new_ct_c2["tissue_ct"].isin(c2_list)]
new_ct_c2 = new_ct_c2.sort_values("Tissue_type")
ct_distance_matrix_c2 = []
for each_ct in new_tissue_c2.Cell_type.unique().tolist():
    print(each_ct)
    ct_matix = new_ct_c2[new_ct_c2["Cell_type"] == each_ct]
    ct_matix = ct_matix.drop(columns=["Cell_type", "tissue_ct"])
    print(ct_matix)
    d_matrix = pairwise_distances(ct_matix.iloc[:, 1:], metric='euclidean')
    print(d_matrix)
    d_list = list(set(i for j in d_matrix for i in j))
    d_list.remove(0)
    print(d_list)
    ct_distance_matrix_c2.append(d_list)
    print()
print([x for xs in ct_distance_matrix_c2 for x in xs])

# Fig.1H
plt.figure().set_figwidth(2.2)

sns.set(style="ticks")
boxplot_data = [[x for xs in ct_distance_matrix_c2 for x in xs], [x for xs in tissue_distance_matrix_c2 for x in xs]]
_ = plt.boxplot(boxplot_data, patch_artist = True, boxprops = dict(facecolor = "green"), medianprops = dict(color = "red"),
                widths=(0.4, 0.4))
_ = plt.xticks([1, 2], ['Tissue type', 'EC type'])
_ = plt.ylabel("Distance")
_ = plt.title("cluster 2")
plt.ylim(0, 27)
plt.tight_layout()
# plt.savefig()
