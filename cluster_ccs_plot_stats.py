# cluster ccs plot and stats

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


def plot_histogram(array, dirpath, label):
    plt.hist(array, bins='auto', color='black', alpha=0.7, label=label)
    plt.xlabel('CCS')
    plt.ylabel('Count')
    plt.savefig(os.path.join(dirpath, 'ccs_distribution_' + label + '.png'))
    plt.close()


def plot_ccs_vs_rosscore(cluster_id, ccs_clust, ros_clust, dirpath):
    color_map = cm.rainbow(np.linspace(0, 1, len(cluster_id)))

    for ind, (clust_id, ccs, ros_score, color) in enumerate(zip(cluster_id, ccs_clust, ros_clust, color_map)):
        plt.scatter(ccs, ros_score, color=color, alpha=0.3, label='cluster '+str(clust_id))

    plt.xlabel('CCS')
    plt.ylabel('Rosetta_score')
    plt.legend(loc='best', fontsize='small')
    plt.savefig(os.path.join(dirpath, 'cluster_ccs_vs_rosscore_distribution.png'), dpi=500)
    plt.close()



def cluster_ccs_stats(dataframe, dirpath, min_clust_size_percent=5.0):

    cluster_id = np.unique(dataframe['clust_num'].values)
    total_size = len(dataframe['clust_num'].values)

    ccs_tm_df = dataframe['ccs_tm'].values
    plot_histogram(ccs_tm_df, dirpath, label='all')

    min_clusters_ccs = []
    min_clusters_rosscore = []
    min_clusters_id = []

    with open(os.path.join(dirpath, 'clust_ccs_stats.csv'), 'w') as outfile:
        header='clust_id,num_structs,percent_structs,mean_ccs,std_ccs\n'
        outfile.write(header)
        for ind, clust_id in enumerate(cluster_id):
            clust_df = dataframe[dataframe['clust_num'] == clust_id]
            clust_df_ccs = clust_df['ccs_tm'].values
            clust_ros_score = clust_df['rosetta_score'].values
            plot_histogram(clust_df_ccs, dirpath, label=str(clust_id))
            mean_ccs = np.mean(clust_df_ccs)
            std_ccs = np.std(clust_df_ccs)
            clust_size = len(clust_df_ccs)
            clust_size_percent = clust_size*100/total_size
            line = '{},{},{},{},{}\n'.format(clust_id, clust_size, clust_size_percent, mean_ccs, std_ccs)
            outfile.write(line)
            if clust_size_percent >= min_clust_size_percent:
                min_clusters_id.append(clust_id)
                min_clusters_ccs.append(clust_df_ccs)
                min_clusters_rosscore.append(clust_ros_score)

        outfile.close()

    plot_ccs_vs_rosscore(min_clusters_id, min_clusters_ccs, min_clusters_rosscore, dirpath)



if __name__=='__main__':
    cluster_pdb_ccs_fpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\cluster_1\num_clusters_6\cluster_pdb_ccs.csv"
    dirpath = os.path.split(cluster_pdb_ccs_fpath)[0]
    df = pd.read_csv(cluster_pdb_ccs_fpath, sep=',')
    cluster_ccs_stats(df, dirpath)
