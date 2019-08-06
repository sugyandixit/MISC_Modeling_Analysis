# select the top k number of interaction from z score combo
#
import os
import pandas as pd
import numpy as np


def get_files(dirpath, endid, startid=None):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith(endid)]
    if startid:
        files = [x for x in files if x.startswith(startid)]
    return files


def combine_top_k_entries(list_of_df, k_entries=10, min_score=8):

    list_df_combine = []

    for ind, df_ in enumerate(list_of_df):

        df_sort = df_.sort_values('comb_score', ascending=False)
        df_new = df_sort[:k_entries]
        df_new = df_new[df_new['comb_score'] > min_score]
        list_df_combine.append(df_new)

    df_combine = pd.concat(list_df_combine)


    return df_combine



def get_combine_df(main_clust_id, dirpath, endid, startid, k_entries=10):

    files = get_files(dirpath, endid, startid)

    list_df = []

    for index, file in enumerate(files):

        fname = str(file).split('.csv')[0]
        clust_ = int(fname.split('_')[-1])

        df = pd.read_csv(os.path.join(dirpath, file))
        df['main_clust_id'] = main_clust_id
        df['sub_clust_id'] = clust_

        list_df.append(df)

    df_combine = combine_top_k_entries(list_df, k_entries=k_entries)

    outname = 'top_k_pair_zscore_combo_all_clust.csv'

    df_combine.to_csv(os.path.join(dirpath, outname), index=False)


def combine_top_k_entries_csv_files(dirpath, endid, startid):

    files = get_files(dirpath, endid, startid)

    list_df = []

    for ind, file in enumerate(files):

        df = pd.read_csv(os.path.join(dirpath, file))
        list_df.append(df)

    df_combine = pd.concat(list_df)

    outname = 'all_clust_top_k_pair_zscore_combo.csv'

    df_combine.to_csv(os.path.join(dirpath, outname), index=False)




if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\top_k_entries_clusters_5"
    endid = '.csv'
    startid = 'top_k_pair_zscore_combo_'
    # get_combine_df(5, dirpath, endid, startid, k_entries=5)
    combine_top_k_entries_csv_files(dirpath, endid, startid)
