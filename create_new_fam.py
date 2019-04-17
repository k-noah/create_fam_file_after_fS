import numpy as np
import pandas as pd
import os
import glob
import sys


'''
use: python create_new_fam.py folder number_of_clusters
'''

folder=sys.argv[1]
sub_fold=folder.split('S')[-1]
no_clust=sys.argv[2]


if int(sub_fold)>1:
    df_fam=pd.read_table('2_mut_miss_'+sub_fold+'.fam',sep=' ',header=None)
    df_clst=pd.read_table(folder+'/2_mut_miss_'+sub_fold+'_output.'+no_clust+'.meanQ',header=None,sep='  ')
else:
    df_fam=pd.read_table('2_mut_miss.fam',sep='\t',header=None)
    df_clst=pd.read_table(folder+'/2_mut_miss_output.'+no_clust+'.meanQ',header=None,sep='  ')
df_clst['Ind']=df_fam[1]

df_clst['max']=''
df_clst['max_clust']=''
for i, row in df_clst.iterrows():
    mx=max(row[:int(no_clust)])
    clst_mx=row[row==mx].index[0]
    df_clst.set_value(i,'max',mx)
    #print(mx)
    if mx>0.66:
        df_clst.set_value(i,'max_clust',clst_mx)

df_cl_pl=pd.DataFrame(columns=['Fam','Ind','Clust'])

for i, row in df_clst.iterrows():
    df_cl_pl.set_value(i,'Ind',row['Ind'])
    df_cl_pl.set_value(i,'Clust',row['max_clust'])

df_cl_pl.replace({0:int(no_clust)},inplace=True)
df_cl_pl['Fam']=0
df_cl_pl.replace('','NA',inplace=True)
df_cl_pl.to_csv(folder+'/2_mut_miss_'+sub_fold+'_cluster.txt',header=False,index=False,sep='\t')