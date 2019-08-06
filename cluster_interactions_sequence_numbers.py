# generate cluster via distance matrix

import os
from math import sqrt
import itertools
import matplotlib.pyplot as plt
import pandas as pd
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as hclust
from scipy.cluster.hierarchy import cut_tree
import numpy as np
from f_test_stats import f_test
import pickle
from scipy.stats import linregress
from plot_landscape_3d import plot_landscape
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
from sklearn.manifold import TSNE


class ClusterDataObject(object):
    """
    save cluster specific pdb file names and pairwise rmsd list
    """
    def __init__(self):
        self.cluster_pdb_files = []
        self.cluster_rmsd_list = []

def save_object_to_pickle(obj, file_path):
    with open(file_path, 'wb') as out_obj:
        pickle.dump(obj, out_obj)
        out_obj.close()


def euclidean_distance_two_dimensions(point1, point2):
    distance = sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2)
    return  distance


def get_distance_list(list_of_points):

    distance_list = []

    for index, (pair1, pair2) in enumerate(itertools.combinations(list_of_points, 2)):

        distance = euclidean_distance_two_dimensions(pair1, pair2)
        distance_list.append(distance)

    return distance_list


def write_anova_ftest_output(ftest, cluster_range, dirpath):
    outlines = ''
    header = '#num_clusters,between_group_var,within_group_var,f_val,p_val\n'
    outlines += header
    for ind in range(len(cluster_range)):
        f_out = ','.join([str(x) for x in ftest[ind]])
        line = '{},{}\n'.format(cluster_range[ind], f_out)
        outlines += line

    with open(os.path.join(dirpath, 'cluster_ftest_output.csv'), 'w') as foutfile:
        foutfile.write(outlines)
        foutfile.close()


def plot_anova_ftest_output(f_val, cluster_range, dirpath):
    plt.scatter(cluster_range, f_val, color='black')
    plt.savefig(os.path.join(dirpath, 'cluster_ftest_output.png'))
    plt.xlabel('Num of clusters')
    plt.ylabel('f-value')
    plt.close()



def f_test_cluster_numbers(linkage, rmsd_sq_mat, dirpath, max_num_clusters=10):
    min_num_clusters=2
    cluster_range = np.arange(min_num_clusters, max_num_clusters+1, 1)
    leaves = hclust.leaves_list(linkage)
    index_list = np.arange(0, len(leaves))
    anova_ftest = []
    test_clusters = []
    test_rmsd_clusts = []
    for num_clust in cluster_range:
        cuttree = cut_tree(linkage, n_clusters=num_clust)
        cuttree_ordered = np.reshape(cuttree, len(leaves))

        clusters = [[] for i in range(num_clust)]
        rmsd_clust = [[] for i in range(num_clust)]

        for ind, leaf in enumerate(leaves):
            cluster = cuttree_ordered[ind]
            ind_num = index_list[ind]
            clusters[cluster].append(ind_num)

        # clusters_ = np.array([np.array(x) for x in clusters])
        test_clusters.append(clusters)

        for num in range(num_clust):
            for index, ind_pair in enumerate(itertools.combinations(clusters[num], 2)):
                rmsd_val = rmsd_sq_mat[ind_pair[0], ind_pair[1]]
                rmsd_clust[num].append(rmsd_val)

        rmsd_clust_ = np.array([np.array(x) for x in rmsd_clust])

        test_rmsd_clusts.append(rmsd_clust)

        between_group_var, within_group_var, f_val, pdf_ = f_test(rmsd_clust_)
        anova_ftest.append([between_group_var, within_group_var, f_val, pdf_])

    anova_ftest = np.array(anova_ftest)

    write_anova_ftest_output(anova_ftest, cluster_range, dirpath)
    plot_anova_ftest_output(anova_ftest[:, 2], cluster_range, dirpath)

    max_fval = np.max(anova_ftest[:, 2])
    ind_max_fval = np.where(anova_ftest[:, 2] == max_fval)[0]
    for ind, (num, test_clust, test_rmsd_clust) in enumerate(zip(cluster_range, test_clusters, test_rmsd_clusts)):
        if ind == ind_max_fval:
            final_num_cluster = num
            final_clusters = test_clust
            final_rmsd_clust = test_rmsd_clust

    return final_num_cluster, final_clusters, final_rmsd_clust



def cluster_distance_list(distance_sq_form):
    linkage = hclust.linkage(distance_sq_form, method='average')
    return linkage


def generate_cluster_data(clusters, cluster_rmsd_list, list_of_pair_points, list_of_pair_residues, list_of_score, dirpath, cluster_mode='auto'):
    clust_data_obj = ClusterDataObject()
    clust_res_pair_store = [[] for i in range(len(clusters))]
    with open(os.path.join(dirpath, 'cluster_' + cluster_mode + '_pair_file.csv'), 'w') as outfile:
        header = 'clust_num,res1_num,res2_num,res1_,res2_,combo_score\n'
        outfile.write(header)
        for index, clust in enumerate(clusters):
            for item in clust:
                pair_points = list_of_pair_points[item]
                pair_point_string = ','.join([str(x) for x in pair_points])
                pair_residues = list_of_pair_residues[item]
                pair_residues_string = ','.join([x for x in pair_residues])
                score_val = list_of_score[item]
                # for ind, pdb_fname in enumerate(rmsd_mat_obj.input_pdb_files):
                #     if index_clust == ind:
                line = '{},{},{},{}\n'.format(index, pair_point_string, pair_residues_string, score_val)
                outfile.write(line)
                clust_res_pair_store[index].append(pair_points)
        outfile.close()
    clust_data_obj.cluster_pdb_files = clust_res_pair_store
    clust_data_obj.cluster_rmsd_list = cluster_rmsd_list
    save_object_to_pickle(clust_data_obj, os.path.join(dirpath, 'cluster_' + cluster_mode + '_data_obj.obj'))


def write_dist_to_rmsd(ht, rmsd_avg_all, rmsd_avg_between_cluster, dirpath, cluster_mode='auto'):
    outlines = ''
    header = '#dist_ht,rmsd_avg_all,rmsd_avg_between_cluster\n'
    outlines += header

    for ind, (height, rmsd_avg_all_, rmsd_avg_between_cluster_) in enumerate(zip(ht, rmsd_avg_all, rmsd_avg_between_cluster)):
        line = '{},{},{}\n'.format(height, rmsd_avg_all_, rmsd_avg_between_cluster_)
        outlines += line

    with open(os.path.join(dirpath, 'cluster_' + cluster_mode + '_dist_to_rmsd.csv'), 'w') as outfile:
        outfile.write(outlines)
        outfile.close()



def dist_to_rmsd(linkage, rmsd_sq_mat, rmsd_dist_mat, dirpath, cluster_mode='auto'):
    dist_range = np.linspace(np.min(rmsd_dist_mat), np.max(rmsd_dist_mat), num=25)
    dist_range = dist_range[1:-2]
    leaves = hclust.leaves_list(linkage)
    index_list = np.arange(0, len(leaves))

    test_rmsd_clusts = []

    for dist in dist_range:
        cuttree = cut_tree(linkage, height=dist)
        cuttree_ordered = np.reshape(cuttree, len(leaves))
        num_clust = np.max(cuttree_ordered) + 1

        clusters = [[] for i in range(num_clust)]
        rmsd_clust = [[] for i in range(num_clust)]

        for ind, leaf in enumerate(leaves):
            cluster = cuttree_ordered[ind]
            ind_num = index_list[ind]
            clusters[cluster].append(ind_num)

        # test_clusters.append(clusters)

        for num in range(num_clust):
            for index, ind_pair in enumerate(itertools.combinations(clusters[num], 2)):
                rmsd_val = rmsd_sq_mat[ind_pair[0], ind_pair[1]]
                rmsd_clust[num].append(rmsd_val)

        test_rmsd_clusts.append(rmsd_clust)

    rmsd_at_height_avg_within_group = []
    rmsd_at_height_avg_all_group = []

    for ind, (height, rmsd_clust_) in enumerate(zip(dist_range, test_rmsd_clusts)):
        # rmsd_clust_ = [np.array(x) for x in rmsd_clust_]
        rmsd_all = np.concatenate(rmsd_clust_)
        mean_rmsd_all = np.nanmean(rmsd_all)
        rmsd_at_height_avg_all_group.append(mean_rmsd_all)

        mean_rmsd_each_group_at_ht = []
        for rmsd_arr in rmsd_clust_:
            try:
                mean_rmsd_each_group_at_ht.append(np.nanmean(rmsd_arr))
            except ValueError:
                pass
        mean_rmsd_each_group_at_ht = np.array(mean_rmsd_each_group_at_ht)
        # mean_rmsd_each_group_at_ht = np.mean(rmsd_clust_, axis=1)
        # rmsd_each_group_at_ht_concatenate = np.concatenate(mean_rmsd_each_group_at_ht)
        mean_rmsd_group_at_ht = np.nanmean(mean_rmsd_each_group_at_ht)
        rmsd_at_height_avg_within_group.append(mean_rmsd_group_at_ht)

    write_dist_to_rmsd(dist_range, rmsd_at_height_avg_all_group, rmsd_at_height_avg_within_group, dirpath,
                       cluster_mode=cluster_mode)

    return dist_range, rmsd_at_height_avg_all_group, rmsd_at_height_avg_within_group

def cut_ht_for_optimal_clusters(dist_range, rmsd_avg_arr, rmsd_avg_opt, upper_lim=-2):
    dist_range = dist_range[:upper_lim]
    rmsd_avg_arr = rmsd_avg_arr[:upper_lim]
    linreg = linregress(rmsd_avg_arr, dist_range)
    cut_height = linreg[0]*rmsd_avg_opt + linreg[1]
    return cut_height




def manual_cluster(dirpath, file, num_clusters=2, cluster_mode='manual', cut_height=None):

    # distance_list, distance_sq_form, list_of_points, list_of_residues, linkage

    df = pd.read_csv(os.path.join(dirpath, file))

    list_of_combo_score = df['comb_score'].values

    list_of_points = []
    list_of_residues = []

    for ind, (res1_num, res2_num, res1, res2) in enumerate(
            zip(df['res1_num'].values, df['res2_num'].values, df['res1_'].values, df['res2_'].values)):
        point = [res1_num, res2_num]
        residues_ = [res1, res2]

        list_of_points.append(point)
        list_of_residues.append(residues_)

    distance_list = get_distance_list(list_of_points)
    distance_sq_form = squareform(distance_list)


    tsne = TSNE(n_components=3, learning_rate=100, metric='precomputed')
    x_tsne = tsne.fit_transform(distance_sq_form)

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    #
    # ax.scatter(x_tsne[:, 0], x_tsne[:, 1], x_tsne[:, 2])
    #
    # plt.show()
    # plt.close()

    print('heho')



    pca_ = PCA(n_components=3)
    pca_.fit(distance_sq_form)

    x_pca = pca_.transform(distance_sq_form)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x_pca[:,0], x_pca[:,1], x_pca[:,2])

    plt.show()
    plt.close()

    print('heho')

    # plot_landscape(distance_sq_form, list_of_combo_score, dirpath, window_savgol=15)

    plt.hist(distance_list, bins='auto', color='grey')
    plt.savefig(os.path.join(dirpath, 'top_k_pair_distance_histogram.png'))
    plt.close()

    linkage = cluster_distance_list(distance_sq_form)


    leaves = hclust.leaves_list(linkage)
    index_list = np.arange(0, len(leaves))

    cuttree = cut_tree(linkage, n_clusters=num_clusters)
    cuttree_ordered = np.reshape(cuttree, len(leaves))

    clusters = [[] for i in range(num_clusters)]
    rmsd_clust = [[] for i in range(num_clusters)]


    for ind, leaf in enumerate(leaves):
        cluster = cuttree_ordered[ind]
        ind_num = index_list[ind]
        clusters[cluster].append(ind_num)


    for num in range(num_clusters):
        for index, ind_pair in enumerate(itertools.combinations(clusters[num], 2)):
            rmsd_val = distance_sq_form[ind_pair[0], ind_pair[1]]
            rmsd_clust[num].append(rmsd_val)

    generate_cluster_data(clusters, rmsd_clust, list_of_points, list_of_residues, list_of_combo_score, dirpath, cluster_mode=cluster_mode)
    cluster_rmsd_avg = np.mean(np.concatenate(rmsd_clust))
    print('Cluster rmsd average: ', cluster_rmsd_avg)
    dist_range, rmsd_avg_all, rmsd_avg_between_clust = dist_to_rmsd(linkage, distance_sq_form, distance_list, dirpath,
                                                                    cluster_mode=cluster_mode)

    if cut_height == None:
        cut_height = cut_ht_for_optimal_clusters(dist_range, rmsd_avg_all, cluster_rmsd_avg, upper_lim=-7)
        print('Cut height: ', cut_height)
    else:
        cut_height = cut_height

    # plot figure

    plt.figure()
    dn = hclust.dendrogram(linkage, distance_sort='ascending', truncate_mode='lastp', color_threshold=cut_height)
    plt.savefig(os.path.join(dirpath, 'top_k_pair_clust_' + cluster_mode + '_dendrogram.png'))
    plt.close()




def get_cluster_from_file_auto(dirpath, file):

    df = pd.read_csv(os.path.join(dirpath, file))

    list_of_points = []
    list_of_residues = []

    list_of_combo_score = df['comb_score'].values


    for ind, (res1_num, res2_num, res1, res2) in enumerate(
            zip(df['res1_num'].values, df['res2_num'].values, df['res1_'].values, df['res2_'].values)):


        point = [res1_num, res2_num]
        residues_ = [res1, res2]

        list_of_points.append(point)
        list_of_residues.append(residues_)

    distance_list = get_distance_list(list_of_points)
    distance_sq_form = squareform(distance_list)

    plot_landscape(distance_sq_form, list_of_combo_score, dirpath)

    plt.hist(distance_list, bins='auto', color='grey')
    plt.savefig(os.path.join(dirpath, 'top_k_pair_distance_histogram.png'))
    plt.close()



    linkage = cluster_distance_list(distance_sq_form)



    plt.figure()
    dn = hclust.dendrogram(linkage, distance_sort='ascending', truncate_mode='lastp')
    plt.savefig(os.path.join(dirpath, 'top_k_pair_clust_auto_dendrogram.png'))
    plt.close()

    num_clusters, clusters, cluster_rmsd = f_test_cluster_numbers(linkage, distance_sq_form, dirpath, max_num_clusters=10)

    generate_cluster_data(clusters, cluster_rmsd, list_of_points, list_of_residues, list_of_combo_score, dirpath)

    print('heho')




if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\top_k_entries_clusters_10"
    fname = "all_clust_top_k_pair_zscore_combo.csv"
    # get_cluster_from_file_auto(dirpath, file=fname)
    manual_cluster(dirpath, fname, num_clusters=4)