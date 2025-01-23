import os
import pandas as pd
import seaborn as sns
import matplotlib.pylab as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# load results from LDA
df = pd.read_csv("./data/lda/mebocost/mebocost_num_norm_receiver_LDA_topics_15t.csv")

important_ms_df = pd.DataFrame()
for i in range(1,16):
    df_topic = df.loc[df["topic"] == i]
    df_topic = df_topic.sort_values("beta", ascending=False)[:10]
    pivot_table = pd.DataFrame(pd.pivot_table(df_topic, values="beta", index="term"))
    pivot_table = pivot_table.rename(columns={"beta": i})
    important_ms_df = pd.concat([important_ms_df, pivot_table], axis=1)

col_order = [1, 11, 5, 9, 3, 7, 2, 8, 13, 12, 4, 6, 15, 10, 14] # based on clustering results in heatmap Fig.4A
row_order = important_ms_df.index.values.tolist()

important_ms_df = important_ms_df.loc[:, col_order]
important_ms_df = important_ms_df.fillna(0)
important_ms_df_norm = important_ms_df.apply(lambda x: (x-x.mean())/x.std(), axis = 1)

# Fig.4B
c_plot = sns.clustermap(important_ms_df_norm, linewidths=0.5, cmap="RdBu_r", col_cluster=False, figsize=(8, 13))
# plt.savefig()

row_order = c_plot.dendrogram_row.reordered_ind
ordered_rows = important_ms_df_norm.index[row_order]

# load in communication results from MEBOCOST
mydata = dict()
results_folder = "./data/mebocost/csv_file"

for filename in sorted(os.listdir(results_folder)):
    filename = os.fsdecode(filename)
    path = results_folder + "/" + filename
    tissue_name = filename.split('_')[0]
    df = pd.read_csv(path, sep='\t')
    mydata[tissue_name] = pd.read_csv(path, sep='\t')
    mydata[tissue_name]["Met-Sensor"] = mydata[tissue_name]["Metabolite_Name"] + "^" + mydata[tissue_name]["Sensor"]

ori = pd.DataFrame(columns=['tissue', 'Met-Sensor', 'Commu_Score'])
for key, df in mydata.items():
    df = df[df['Receiver'].str.contains('EC') == True]
    num = pd.DataFrame(df.groupby('Met-Sensor').count())
    num['tissue'] = key
    num.reset_index(inplace=True)
    num = num.rename(columns = {'index':'Met-Sensor'})
    ori = pd.concat([ori, num], axis=0)
    print(ori)
new = ori.copy()
ori = ori.pivot(columns='tissue', values='Commu_Score', index='Met-Sensor')
new = new.pivot(columns='tissue', values='Commu_Score', index='Met-Sensor')
new = new.fillna(0)

new2 = new[new.index.isin(ordered_rows)]
new2 = new2.loc[ordered_rows, :]
new2_norm= new2.apply(lambda x: (x-x.mean())/x.std(), axis = 1)
col_order2 = ["testis", "thymus", "bladder", "liver", "skin", "ovary", "adipose", "breast", "heart", "kidney", "skeletalMuscle", "stomach", "lung", "gut", "intestine"] # order based on Fig.4A
new = new.loc[row_order, col_order2]

# Fig.4C
ax, fig = plt.subplots(figsize=(7,9))
sns.heatmap(new2_norm, linewidths=0.5, cmap="RdBu_r", cbar_kws={'label': 'normalized communication number'})
plt.xlabel("Tissue")
plt.ylabel("Metabolite-Sensor Pair")
plt.title("MEBOCOST")
plt.tight_layout()
# plt.savefig()

ms_interest_filtered_mydata = mydata.copy()
# select ms pairs of interests
met_sensors_of_interest = ["Cholesterol^RORA", "L-Lysine^SLC7A1", "Ornithine^SLC7A1", "Adenine^SLC35F5", "L-Cysteine^SLC43A2"]
for tissue, df in ms_interest_filtered_mydata.items():
    ms_interest_filtered_mydata[tissue] = df.loc[df["Met-Sensor"].isin(met_sensors_of_interest)]
combined_data = pd.concat([df.assign(Tissue=tissue) for tissue, df in ms_interest_filtered_mydata.items()])

unique_met_sensors = met_sensors_of_interest
n_met_sensors = len(unique_met_sensors)
n_tissues = len(mydata)
combined_data["Tissue"] = pd.Categorical(combined_data["Tissue"], categories=sorted(combined_data["Tissue"].unique()), ordered=True)

# Fig.4D
fig, axes = plt.subplots(n_met_sensors, 1, figsize=(10, 9), sharex=True)
for ax, met_sensor in zip(axes, unique_met_sensors):
    norm = mcolors.Normalize(vmin=combined_data["Commu_Score"].min(), vmax=combined_data["Commu_Score"].max())
    colormap = cm.get_cmap('Greens')
    colors = colormap(norm(combined_data["Commu_Score"].values))
    tissue_colors = dict(zip(combined_data["Tissue"], colors)) 
    sns.violinplot(
        x="Tissue", y="Commu_Score", data=combined_data[combined_data["Met-Sensor"] == met_sensor],
        ax=ax, inner=None, palette=tissue_colors
    )
    ax.set_title(f"Met-Sensor: {met_sensor}", fontsize=14)
    ax.set_xlabel("Tissue", fontsize=14)
    ax.set_ylabel("Communication Score", fontsize=14)
sm = plt.cm.ScalarMappable(cmap=colormap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes, orientation='vertical', fraction=0.02, pad=0.02, location="right")
cbar.set_label('Communication Score', fontsize=14)
# plt.savefig()

