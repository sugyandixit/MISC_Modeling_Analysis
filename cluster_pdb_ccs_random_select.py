# randomly select n number of structures from each cluster and make an output file
# can be used to select a subset of structures for trajectory method ccs calculation

import os
import pandas as pd
import numpy as np

def sel_randomly_from_each_cluster(filepath, sample_num=5):
    """
    select randomly from each cluster given the number to sample
    :param filepath: file path to choose from
    :param sample_num: number to select randomly
    :return: void. Outputfile
    """

    df = pd.read_csv(filepath)
    uniq_clusts = np.unique(df['clust_num'].values)

    header = 'clust_num,pdb_fname,CCS_PA,CCS_TM\n'
    data_string = ''


    for ind, clust in enumerate(uniq_clusts):
        df_clust = df[df['clust_num'] == clust]
        len_clust = len(df_clust['clust_num'].values)
        if len_clust > sample_num:
            random_ind = np.random.choice(len_clust-1, sample_num)
            for random_index in random_ind:
                df_rand = df_clust.iloc[[random_index]]
                line = '{},{},{},{}\n'.format(df_rand['clust_num'].values[0],
                                              df_rand['pdb_fname'].values[0],
                                              df_rand['CCS_PA'].values[0],
                                              df_rand['CCS_TM'].values[0])
                data_string += line
        else:
            for num, (clust_num, pdb_fname, ccs_pa, ccs_tm) in enumerate(zip(df_clust['clust_num'].values,
                                                                             df_clust['pdb_fname'].values,
                                                                             df_clust['CCS_PA'].values,
                                                                             df_clust['CCS_TM'].values)):
                line = '{},{},{},{}\n'.format(clust_num,
                                              pdb_fname,
                                              ccs_pa,
                                              ccs_tm)
                data_string += line

    output_string = header + data_string

    dirpath, orig_fname = os.path.split(filepath)
    out_fname = str(orig_fname).split('.csv')[0] + '_randomselection.csv'

    with open(os.path.join(dirpath, out_fname), 'w') as outfile:
        outfile.write(output_string)
        outfile.close()

    print('heho')



if __name__=='__main__':

    clust_pdb_ccs_fpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_nomod\rmsdmat\cluster_manual_pdb_ccs.csv"
    sel_randomly_from_each_cluster(clust_pdb_ccs_fpath, sample_num=10)