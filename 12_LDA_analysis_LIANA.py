import os
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
from statsmodels.stats.multitest import fdrcorrection
from itertools import product

# load results from LDA
df = pd.read_csv("./data/lda/liana/liana_num_norm_receiver_LDA_topics_15t.csv")

important_lr_df = pd.DataFrame()
for i in range(1,16):
    df_topic = df.loc[df["topic"] == i]
    df_topic = df_topic.sort_values("beta", ascending=False)[:10]
    pivot_table = pd.DataFrame(pd.pivot_table(df_topic, values="beta", index="term"))
    pivot_table = pivot_table.rename(columns={"beta": i})
    important_lr_df = pd.concat([important_lr_df, pivot_table], axis=1)
important_lr_df = important_lr_df.fillna(0)

col_order = [1, 7, 13, 15, 2, 10, 11, 12, 8, 4, 14, 5, 6, 3, 9] # order based on clustering results in heatmap Fig.5A
row_order = important_lr_df.index.values.tolist()
important_lr_df = important_lr_df.loc[:, col_order]
important_lr_df_norm = important_lr_df.apply(lambda x: (x-x.mean())/x.std(), axis = 1)

# Fig.5B
c_plot = sns.clustermap(important_lr_df_norm, linewidths=0.5, cmap="RdBu_r", col_cluster=False, figsize=(8, 19))
# plt.savefig()

row_order = c_plot.dendrogram_row.reordered_ind
ordered_rows = important_lr_df_norm.index[row_order]

# load in communication results from LIANA
mydata = dict()
results_folder = "./data/liana"

for filename in sorted(os.listdir(results_folder)):
    filename = os.fsdecode(filename)
    path = results_folder + "/" + filename
    tissue_name = filename.split('_')[0]
    mydata[tissue_name] = pd.read_csv(path, sep=',')
    fdrcorrection_res = fdrcorrection(mydata[tissue_name]["specificity_rank"].to_numpy(), alpha=0.05, method='indep', is_sorted=False)
    fdrcorrection_res = pd.Series(fdrcorrection_res)
    mydata[tissue_name].insert(13, "fdrcorrection", fdrcorrection_res[1])
    mydata[tissue_name] = mydata[tissue_name][mydata[tissue_name]["fdrcorrection"] >= 0.05]
    mydata[tissue_name]["Ligand-Receptor"] = mydata[tissue_name]["ligand_complex"] + "^" + mydata[tissue_name]["receptor_complex"]

ori2 = pd.DataFrame(columns=['tissue', 'Ligand-Receptor', 'num'])
for key, df in mydata.items():
    df = df[df['target'].str.contains('EC') == True]
    num = pd.DataFrame(df.groupby('Ligand-Receptor').count())
    num['tissue'] = key
    num["num"] = num["source"]
    num.reset_index(inplace=True)
    num = num.rename(columns = {'index':'Ligand-Receptor'})
    ori2 = pd.concat([ori2, num], axis=0)

new2 = ori2.copy()
ori2 = ori2.pivot(columns='tissue', values='num', index='Ligand-Receptor')
new2 = new2.pivot(columns='tissue', values='num', index='Ligand-Receptor')

new2 = new2.fillna(0)
new2 = new2[new2.index.isin(ordered_rows)]
new2 = new2.loc[ordered_rows, :]
new2_norm= new2.apply(lambda x: (x-x.mean())/x.std(), axis = 1)
col_order2 = ["thymus", "heart", "skeletalMuscle", "gut", "intestine", "kidney", "stomach", "breast", "adipose", "bladder", "ovary", "liver", "testis", "lung", "skin"] # order based on Fig.5A
new2_norm = new2_norm.loc[ordered_rows, col_order2]

# Fig.5C
ax, fig = plt.subplots(figsize=(7,13))
sns.heatmap(new2_norm, linewidths=0.5, cmap="RdBu_r", cbar_kws={'label': 'normalized communication number'})
plt.xlabel("Tissue")
plt.ylabel("Ligand-Receptor Pair")
plt.title("LIANA")
plt.tight_layout()
# plt.savefig()


ori3 = pd.DataFrame(columns=['tissue', 'Ligand-Receptor', 'lr_means'])

all_ligand_receptors = pd.concat([df[df['target'].str.contains('EC')]['Ligand-Receptor'] for df in mydata.values()]).unique()
for key, df in mydata.items():
    df = df[df['target'].str.contains('EC') == True]
    num = pd.DataFrame(df[["Ligand-Receptor", "lr_means"]]).reset_index()
    all_combinations = pd.DataFrame(list(product([key], all_ligand_receptors)), columns=['tissue', 'Ligand-Receptor'])
    num_full = pd.merge(all_combinations, num, on='Ligand-Receptor', how='left')
    num_full['lr_means'] = num_full['lr_means'].fillna(0)
    ori3 = pd.concat([ori3, num_full], axis=0)
ori3.reset_index(drop=True, inplace=True)

# select lr paris of interests
ligand_receptors = ["ADAM17^ITGB1", "CALM1^RYR2", "ADAM17^ITGA5", "LAMA2^RPSA", "ADAM17^NOTCH1", "MAML2^NOTCH1"] 

# Fig.5D
fig, axes = plt.subplots(nrows=6, ncols=1, figsize=(10, 10)) 
for i, lr in enumerate(ligand_receptors):
    lr_pair = ori3[ori3["Ligand-Receptor"] == lr]
    sns.violinplot(ax=axes[i], data=lr_pair, x="tissue", y="lr_means", inner=None, linewidth=2, scale='width', width=0.8)
    axes[i].set_title(f'Ligand-Receptor: {lr}')
    axes[i].set_xlabel('Tissue')
    axes[i].set_ylabel('lr_means')
plt.tight_layout()
# plt.savefig()
