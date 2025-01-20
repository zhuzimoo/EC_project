import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from scipy import stats

# import MEBOCOST results
mydata_mebocost = dict()
results_folder = "./data/mebocost/csv_file"

for filename in sorted(os.listdir(results_folder)):
    filename = os.fsdecode(filename)
    path = results_folder + "/" + filename
    tissue_name = filename.split('_')[0]
    df = pd.read_csv(path, sep='\t')
    mydata_mebocost[tissue_name] = pd.read_csv(path, sep='\t')

new = pd.DataFrame(columns=['tissue', 'num', 'property'])
for key, value in mydata_mebocost.items():
    tissue = key
    value = value[(value['Sender'] == 'Arterial EC') | (value['Sender'] == 'Venous EC') | (value['Sender'] == 'Capillary EC') | (value['Sender'] == 'Lymphatic EC')| 
                  (value['Receiver'] == 'Arterial EC') | (value['Receiver'] == 'Venous EC') | (value['Receiver'] == 'Capillary EC') | (value['Receiver'] == 'Lymphatic EC')]
    num_commu = value.shape[0]
    new = new._append(pd.Series({'tissue': tissue, 'num': num_commu, 'property': 'mebocost'}), ignore_index=True)

# import LIANA results
mydata_liana = dict()
results_folder = "./data/liana"
for filename in sorted(os.listdir(results_folder)):
    filename = os.fsdecode(filename)
    path = results_folder + "/" + filename
    tissue_name = filename.split('_')[0]
    df = pd.read_csv(path, sep=',')
    mydata_liana[tissue_name] = pd.read_csv(path, sep=',')

# make a new, combined dataframe
for key, value in mydata_liana.items():
    tissue = key
    value = value[(value['source'] == 'Arterial EC') | (value['source'] == 'Venous EC') | (value['source'] == 'Capillary EC') | (value['source'] == 'Lymphatic EC')| 
                  (value['target'] == 'Arterial EC') | (value['target'] == 'Venous EC') | (value['target'] == 'Capillary EC') | (value['target'] == 'Lymphatic EC')]
    num_commu = value.shape[0]
    new = new._append(pd.Series({'tissue': tissue, 'num': num_commu, 'property': 'liana'}), ignore_index=True)

total = new.groupby('tissue')['num'].sum()
order_total = total.sort_values()
sorted_total_list = list(order_total.keys())

# info got from adata
ct_count_dict = {'adipose': 29,
                 'bladder': 36,
                 'breast': 37,
                 'gut': 62,
                 'heart': 29,
                 'intestine': 67,
                 'kidney': 18,
                 'liver': 52,
                 'lung': 45,
                 'ovary': 44,
                 'skeletalMuscle': 30,
                 'skin': 38,
                 'stomach': 29,
                 'testis': 40,
                 'thymus': 49}

total_ct_count = []
for row in order_total.keys():
    total_ct_count.append(ct_count_dict[row])
total_avg_ct = sum(total_ct_count) / len(total_ct_count)

pct_dict = {'adipose': 0.11,
            'bladder': 18.7,
            'breast': 9.89,
            'gut': 2.1999999999999997,
            'heart': 13.900000000000002, 
            'intestine': 2.18,
            'kidney': 3.2399999999999998, 
            'liver': 0.54, 'lung': 6.87,
            'ovary': 3.4099999999999997, 
            'skeletalMuscle': 23.02, 
            'skin': 14.77, 
            'stomach': 0.5662514156285391, 
            'testis': 3.9600000000000004, 
            'thymus': 3.02}

total_pct = []
for row in order_total.keys():
    total_pct.append(pct_dict[row])
total_avg_pct = sum(total_pct) / len(total_pct)

com_per_ct = new.copy()
result = []
for index, row in com_per_ct.iterrows():
    result.append(round(row['num'] / ct_count_dict[row['tissue']], 1))
com_per_ct["com_per_ct"] = result

mebocost_new = new[new["property"] == "mebocost"]
liana_new = new[new['property'] == "liana"]
mebocost_com_per_ct = com_per_ct[com_per_ct['property'] == "mebocost"]
liana_com_per_ct = com_per_ct[com_per_ct['property'] == "liana"]

mebocost_new.reindex(sorted_total_list)
liana_new.reindex(sorted_total_list)
mebocost_new.sort_values(by="tissue", key=lambda column: column.map(lambda e: sorted_total_list.index(e)), inplace=True)
liana_new.sort_values(by="tissue", key=lambda column: column.map(lambda e: sorted_total_list.index(e)), inplace=True)

# Fig.3A
fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(111)
ax2 = ax.twinx()
width = 0.4
p1 = mebocost_new.plot(x='tissue', y='num', kind='bar', color='tab:blue', ax=ax, width=width, position=1, label="metabolite-sensor")
p2 = liana_new.plot(x='tissue', y='num', kind='bar', color='tab:orange', ax=ax2, width=width, position=0, label="ligand-receptor")
h1, l1 = p1.get_legend_handles_labels()
h2, l2 = p2.get_legend_handles_labels()
plt.legend(h1+h2, l1+l2, loc=2)
ax.set_ylabel('metabolite-sensor')
ax2.set_ylabel('ligand-receptor')
plt.tight_layout()
# plt.savefig()

# Fig.3B
fig, axis = plt.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios': [1.5, 1.5]}, figsize=(10, 3))
f2 = sns.barplot(x=np.arange(len(total_pct)), y=total_pct, ax=axis[0], color="yellowgreen")
f2.set(xticklabels=[]) 
f2.set(ylabel="EC pct (%)")
plt.tight_layout()
f2.axhline(y = total_avg_pct, color = 'r', linestyle = '--', label=f"average={round(total_avg_pct,2)}%", alpha=0.5)
f2.legend()
axis[0].bar_label(axis[0].containers[0])
axis[0].set_ylim(0, 35)
axis[0].set_ylabel("EC pct (%)", fontsize = 11)

f3 = sns.barplot(x=np.arange(len(total_ct_count)), y=total_ct_count, ax=axis[1], color="gold")
f3.set(xticklabels=[]) 
f3.set(ylabel="# cell type")
f3.axhline(y = total_avg_ct, color = 'r', linestyle = '--', label=f"average={round(total_avg_ct,2)}", alpha=0.5)
f3.legend()
axis[1].bar_label(axis[1].containers[0])
axis[1].set_ylim(0, 95)
axis[1].set_ylabel("# cell type", fontsize = 11)
# plt.savefig()

# Regression
mydf = new.groupby('tissue', as_index = False)['num'].sum()
mydf['newtissue'] = mydf['tissue'].apply(lambda x: sorted_total_list.index(x))
mydf = mydf.sort_values(by='newtissue')
mydf = mydf.drop(columns=['newtissue'])
mydf['pct'] = total_pct
mydf['ct_count'] = total_ct_count
mydf['num'] = mydf['num'].astype(float)
corr1, _ = pearsonr(mydf['num'], mydf['pct'])
corr2, _ = pearsonr(mydf['num'], mydf['ct_count'])

# Fig.3C
slope, intercept, r_value, p_value, std_err = stats.linregress(mydf['num'],mydf['pct'])
plt.figure(figsize=(5, 5))
sns.lmplot(data=mydf, x='num', y='pct', line_kws={'label':f'correlation coefficient = {round(corr1,2)}'})
plt.legend(fontsize=12)
plt.ylabel("EC pct (%)", fontsize=12)
plt.xlabel("# Communications", fontsize=12)
plt.legend(loc='upper right')
# plt.savefig()

# Fig.3D
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(mydf['num'],mydf['ct_count'])
plt.figure(figsize=(5, 5))
sns.lmplot(data=mydf, x='num', y='ct_count', line_kws={'label':f'correlation coefficient = {round(corr2, 2)}'})
plt.legend(fontsize=12)
plt.ylabel("# cell type", fontsize=12)
plt.xlabel("# Communications", fontsize=12)
plt.legend(loc='upper left')
# plt.savefig()

# Fig.3E
mebocost_com_per_ct.reindex(sorted_total_list)
liana_com_per_ct.reindex(sorted_total_list)
mebocost_com_per_ct.sort_values(by="tissue", key=lambda column: column.map(lambda e: sorted_total_list.index(e)), inplace=True)
liana_com_per_ct.sort_values(by="tissue", key=lambda column: column.map(lambda e: sorted_total_list.index(e)), inplace=True)

fig = plt.figure(figsize=(10, 4))
ax = fig.add_subplot(111)
ax2 = ax.twinx()
width = 0.4
p1 = mebocost_com_per_ct.plot(x='tissue', y='com_per_ct', kind='bar', color='tab:blue', ax=ax, width=width, position=1, label="metabolite-sensor")
p2 = liana_com_per_ct.plot(x='tissue', y='com_per_ct', kind='bar', color='tab:orange', ax=ax2, width=width, position=0, label="ligand-receptor")
h1, l1 = p1.get_legend_handles_labels()
h2, l2 = p2.get_legend_handles_labels()
ax.get_legend().remove()
plt.legend(h1+h2, l1+l2, loc=2)
ax.set_ylabel('metabolite-sensor')
ax2.set_ylabel('ligand-receptor')
ax.set_xlabel("Tissue Type", fontsize=13)
ax.tick_params(axis='x', labelrotation = 90, labelsize=14)
plt.tight_layout()
# plt.savefig()