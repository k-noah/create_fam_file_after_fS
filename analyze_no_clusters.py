import numpy as np
import pandas as pd
import os
import glob
import sys

'''
usage: python analyze_no_clusters.py folder
'''

df_genes=pd.read_table('2_mut_miss.ped',sep='\t',header=None,usecols=range(10))
#df_genes.set_value(1,1,'0C2E_isolate_RP1_3623')

df_cf=pd.read_csv('combined_snps_gt1per_rem_250231_rem-1.csv',index_col=0)
df_cfT=df_cf.T
cf_fl=df_cfT['CF_flag']
l_strains=df_genes[1]

folder=sys.argv[1]
#print(folder)

files=glob.glob(folder+'/*.meanQ')
df_compl=pd.DataFrame()
for file in files:
    print(file)
    sht_file=file.split('.')[-2]
    file_ml=glob.glob(folder+'/*.'+sht_file+'.log')[0]
    df=pd.read_table(file,sep='  ',header=None)
    df=df.set_index(l_strains)
    df['CF_flag']=cf_fl
    df['cluster']=''
    df['no_clst']=0
    for i, row in df.iterrows():
        #print(i)
        l_rw=len(row)
        rw=row[:l_rw-3]
        mx=max(rw)
        col=''
        l_col=0
        if mx>0.6:
            col=str(rw.idxmax(axis=1))+' '
            df.set_value(i,'cluster',col)
            l_col=1
        else:
            while mx>0.3:
                col=col+str(rw.idxmax(axis=1))+' '
                rw=rw.drop(rw.idxmax(axis=1))
                mx=max(rw)
                l_col+=1
            df.set_value(i,'cluster',col)
        df.set_value(i,'no_clst',l_col)
        #pivot table of #clusters
    table = pd.pivot_table(df, values=0, index=['no_clst'], aggfunc='count')
    #table.to_csv(file+'_table1.csv')
    #print(table)
    #pivot table of #CF/non-CF per cluster
    table2 = pd.pivot_table(df, values=0, index=['cluster'], columns=['CF_flag'], aggfunc='count')
    #table2.to_csv(file+'_table2.csv')
    #print(len(table2))
    df.to_csv(file+'_overview.csv')
    for j, row in table.iterrows():
        df_compl.set_value(j,sht_file,row[0])
    df_compl.set_value('tot_clst',sht_file,len(table2))
    with open(file_ml,'r') as f:
        lines=f.readlines()
        ml=lines[-3].split(' ')[-1]
    df_compl.set_value('marginal_likelihood',sht_file,ml)
df_compl.to_csv(folder+'/clusters_overview.csv')