import os
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

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

ori = pd.DataFrame(columns=['tissue', 'Ligand-Receptor', 'num'])
for key, df in mydata.items():
    df = df[df['target'].str.contains('EC') == True]
    num = pd.DataFrame(df.groupby('Ligand-Receptor').count())
    num['tissue'] = key
    num["num"] = num["source"]
    num.reset_index(inplace=True)
    num = num.rename(columns = {'index':'Ligand-Receptor'})
    ori = pd.concat([ori, num], axis=0)

new = ori.copy()
ori = ori.pivot(columns='tissue', values='num', index='Ligand-Receptor')
new = new.pivot(columns='tissue', values='num', index='Ligand-Receptor')
new = new.fillna(0)

new['adipose'] = new['adipose'].apply(lambda x: x*100/29)
new['bladder'] = new['bladder'].apply(lambda x: x*100/36)
new['breast'] = new['breast'].apply(lambda x: x*100/37)
new['gut'] = new['gut'].apply(lambda x: x*100/62)
new['heart'] = new['heart'].apply(lambda x: x*100/29)
new['intestine'] = new['intestine'].apply(lambda x: x*100/67)
new['kidney'] = new['kidney'].apply(lambda x: x*100/18)
new['liver'] = new['liver'].apply(lambda x: x*100/52)
new['lung'] = new['lung'].apply(lambda x: x*100/45)
new['ovary'] = new['ovary'].apply(lambda x: x*100/44)
new['skeletalMuscle'] = new['skeletalMuscle'].apply(lambda x: x*100/30)
new['skin'] = new['skin'].apply(lambda x: x*100/38)
new['stomach'] = new['stomach'].apply(lambda x: x*100/29)
new['testis'] = new['testis'].apply(lambda x: x*100/40)
new['thymus'] = new['thymus'].apply(lambda x: x*100/49)
new = new.round(0)

# new.to_csv()