import os
import pandas as pd

mydata = dict()
# load results from MEBOCOST inferred CCC
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