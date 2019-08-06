# compile a file with ccs appended to each pdb file in the cluster pdb file.csv

import os
import pandas as pd
import numpy as np


impact_ccs_fpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\cluster_1\impact_ccs.csv"
ros_score_fpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\cluster_1\pdb_list_ct_map_5.csv"
cluster_pdb_fpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\cluster_1\num_clusters_6\cluster_pdb_files.csv"

dirpath = os.path.split(cluster_pdb_fpath)[0]

impact_ccs_df = pd.read_csv(impact_ccs_fpath, sep=',')
clust_pdb_df = pd.read_csv(cluster_pdb_fpath, sep=',')
clust_pdb_df = clust_pdb_df.drop_duplicates(subset='pdb_fname')
pdb_rosscore_df = pd.read_csv(ros_score_fpath, sep=',')

clust_pdb_ccs_pa = []
clust_pdb_ccs_tm = []
ros_score = []
for index, clust_pdb_fpath in enumerate(clust_pdb_df['pdb_fname'].values):
    print(index)
    clust_pdb_fname = os.path.split(clust_pdb_fpath)[1]
    for ind, (ccs_pdb_name, ccs_pa, ccs_tm) in enumerate(zip(impact_ccs_df['filename'].values,
                                                             impact_ccs_df['CCS_PA'].values,
                                                             impact_ccs_df['CCS_TM'].values)):
        if clust_pdb_fname == ccs_pdb_name:
            clust_pdb_ccs_pa.append(ccs_pa)
            clust_pdb_ccs_tm.append(ccs_tm)
    rosetta_sc = pdb_rosscore_df[pdb_rosscore_df['pdb_name'] == clust_pdb_fname]['rosetta_score'].values
    ros_score.append(rosetta_sc)

clust_pdb_ccs_pa = np.array(clust_pdb_ccs_pa)
clust_pdb_ccs_tm = np.array(clust_pdb_ccs_tm)
ros_score = np.array(ros_score)

clust_pdb_df['ccs_pa'] = clust_pdb_ccs_pa
clust_pdb_df['ccs_tm'] = clust_pdb_ccs_tm
clust_pdb_df['rosetta_score'] = ros_score

outfile = os.path.join(dirpath, 'cluster_pdb_ccs.csv')
clust_pdb_df.to_csv(outfile, sep=',', index=False)

print('heho')