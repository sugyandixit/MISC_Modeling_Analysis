# combine z score by calculating a combined score
# need three (or more) zscores
# score = wij * (wi + wj)

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def get_files(dirpath, endid, startid=None):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith(endid)]
    if startid:
        files = [x for x in files if x.startswith(startid)]
    return files

def calc_comb_z_scores(pair_zscore, ind_z_score):
    sum_weights = 0
    for ind_score in ind_z_score:
        sum_weights += ind_score
    score = pair_zscore * sum_weights
    return score


def comb_z_score(pair_z_score_df, res1_zscore_df, res2_zscore_df):

    comb_score_arr = []

    for ind, (res1_num, res2_num, z_score) in enumerate(
        zip(pair_z_score_df['res1_num'].values, pair_z_score_df['res2_num'].values, pair_z_score_df['zscore'].values)):
        res1_z_score = res1_zscore_df[res1_zscore_df['residue_num'] == res1_num+1]['zscore'].values[0]
        res2_z_score = res2_zscore_df[res2_zscore_df['residue_num'] == res2_num + 1]['zscore'].values[0]
        comb_score = calc_comb_z_scores(z_score, [res1_z_score, res2_z_score])
        comb_score_arr.append(comb_score)

    pair_z_score_df['comb_score'] = np.array(comb_score_arr)

    return pair_z_score_df


def construct_grid(x, y, data):
    x_unique = np.unique(x)
    y_unique = np.unique(y)
    data_grid = np.zeros(shape=(len(x_unique), len(y_unique)))
    ind_new = 0
    for ind, (x_, y_) in enumerate(zip(x, y)):
        for i, x_un in enumerate(x_unique):
            for j, y_un in enumerate(y_unique):
                if x_un == x_:
                    if y_un == y_:
                        data_grid[i, j] = data[ind_new]
                        ind_new += 1

    data_grid[data_grid == 0] = 'nan'
    return data_grid


def plot_zscore_imhow(zscoremat, cmap, outid, dirpath, xaxis, yaxis, xspacing=1, yspacing=1, image_type='.png'):
    y_ticks = np.arange(0, len(yaxis), yspacing)
    y_tick_labels = []
    for ind in y_ticks:
        y_tick_labels.append(yaxis[ind])
    plt.imshow(zscoremat, origin=('lower', 'lower'), cmap=cmap)
    plt.xticks(np.arange(0, len(xaxis), xspacing), xaxis)
    plt.yticks(y_ticks, y_tick_labels)
    plt.colorbar()
    plt.savefig(os.path.join(dirpath, 'comb_z_score_'+ outid +image_type), dpi=500)
    plt.close()



def comb_z_score_for_clust(dirpath, clust_num, pair_id, res1_id, res2_id):
    pair_zscore_file = get_files(dirpath, endid='zscore_clust_'+str(clust_num)+'.csv', startid=pair_id)
    res1_zscore_file = get_files(dirpath, endid='zscore_cluster_'+str(clust_num)+'.csv', startid=res1_id)
    res2_zscore_file = get_files(dirpath, endid='zscore_cluster_'+str(clust_num)+'.csv', startid=res2_id)

    for index, (pair_file, res1_file, res2_file) in enumerate(zip(pair_zscore_file, res1_zscore_file, res2_zscore_file)):

        pair_df = pd.read_csv(os.path.join(dirpath, pair_file))
        res1_df = pd.read_csv(os.path.join(dirpath, res1_file))
        res2_df = pd.read_csv(os.path.join(dirpath, res2_file))

        comb_z_score_df = comb_z_score(pair_df, res1_df, res2_df)
        comb_z_score_df_sorted = comb_z_score_df.sort_values('comb_score', ascending=False)

        grid_score = construct_grid(pair_df['res2_num'].values, pair_df['res1_num'].values,
                                    comb_z_score_df['comb_score'].values)
        plot_zscore_imhow(grid_score, cmap='PuBu', outid='clust_'+ str(clust_num), dirpath=dirpath,
                          xaxis=res1_df['residue_num'].values,
                          yaxis=res2_df['residue_num'].values)
        comb_z_score_df_sorted.to_csv(os.path.join(dirpath, 'pair_zscore_combo_clust_'+str(clust_num)+'.csv'), index=False)

        # print('hiho')



if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\cluster_5"

    clust_nums = np.arange(0, 5)

    for ind, clust_ in enumerate(clust_nums):

        comb_z_score_for_clust(dirpath, clust_num=clust_, pair_id='pair', res1_id='mol1', res2_id='mol2')