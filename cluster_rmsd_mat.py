# cluster_rmsd_mat.py

import os
import pickle
import itertools
import f_test_stats
from f_test_stats import f_test
import rmsd_mat_obj
from rmsd_mat_obj import RMSDMAT
import numpy as np
from scipy.spatial.distance import squareform, pdist
import matplotlib.pyplot as plt
from matplotlib.pyplot import  cm
from scipy.cluster.hierarchy import cophenet, cut_tree
import scipy.cluster.hierarchy as hclust
from scipy.stats import linregress


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



def generate_cluster_data(clusters, cluster_rmsd_list, rmsd_mat_obj, dirpath, cluster_mode='auto'):
    clust_data_obj = ClusterDataObject()
    clust_pdb_store = [[] for i in range(len(clusters))]
    with open(os.path.join(dirpath, 'cluster_' + cluster_mode + '_pdb_files.csv'), 'w') as outfile:
        header = 'clust_num,pdb_fname\n'
        outfile.write(header)
        for index, clust in enumerate(clusters):
            for item in clust:
                pdb_fname = rmsd_mat_obj.input_pdb_files[item]
                # for ind, pdb_fname in enumerate(rmsd_mat_obj.input_pdb_files):
                #     if index_clust == ind:
                line = '{},{}\n'.format(index, pdb_fname)
                outfile.write(line)
                clust_pdb_store[index].append(pdb_fname)
        outfile.close()
    clust_data_obj.cluster_pdb_files = clust_pdb_store
    clust_data_obj.cluster_rmsd_list = cluster_rmsd_list
    save_object_to_pickle(clust_data_obj, os.path.join(dirpath, 'cluster_' + cluster_mode + '_data_obj.obj'))
    # return


def load_pickle_object(pickle_fpath):
    """
        opens the pickled object!
        Caution: You have to import the class(es) into current python script in order to load the object
        First import the python script: import python_script
        then import the class: from python_script import class(es)
        :param pickle_file_path: pickle object path
        :return: object
        """
    with open(pickle_fpath, 'rb') as pk_file:
        obj = pickle.load(pk_file)
    return obj


def plot_rmsd_dist_mat_scatter(rmsd_mat_sqform, dirpath):
    plt.scatter(rmsd_mat_sqform[:, 0], rmsd_mat_sqform[:, 1])
    plt.savefig(os.path.join(dirpath, rmsd_mat_obj_fname + '_rmsd_dist_scatter_plot.png'))
    plt.close()

def plot_rmsd_dist_hist(rmsd_mat, dirpath):
    plt.hist(rmsd_mat, bins='auto', color='black', alpha=0.7)
    plt.savefig(os.path.join(dirpath, rmsd_mat_obj_fname + '_rmsd_dist_hist.png'))
    plt.close()


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


# def generate_color_link_palette_with_num_clusters(num_clusters, cmap):
#     """
#     provide the number of clusters and return the
#     :param num_clusters: number of clusters
#     :param cmap: color map
#     :return:link color list
#     """
#     uncluster_col = '#808080'
#
#     d_leaf_



def autogen_cluster(rmsd_mat_obj, rmsd_sq_mat, linkage, dirpath, max_num_clusters = 10, cluster_mode='auto', cut_height=None):
    num_clusters, clusters, cluster_rmsd = f_test_cluster_numbers(linkage, rmsd_sq_mat, dirpath, max_num_clusters)
    generate_cluster_data(clusters, cluster_rmsd, rmsd_mat_obj, dirpath, cluster_mode=cluster_mode)
    cluster_rmsd_avg = np.mean(np.concatenate(cluster_rmsd))
    print('Cluster rmsd average: ', cluster_rmsd_avg)
    dist_range, rmsd_avg_all, rmsd_avg_between_clust = dist_to_rmsd(linkage, rmsd_mat_sqform, rmsd_dist_mat,
                                                                    dirpath, cluster_mode=cluster_mode)

    if cut_height == None:
        cut_height = cut_ht_for_optimal_clusters(dist_range, rmsd_avg_all, cluster_rmsd_avg)
        print('Cut height: ', cut_height)
    else:
        cut_height = cut_height

    #plot figure

    plt.figure()
    dn = hclust.dendrogram(linkage, distance_sort='ascending', truncate_mode='lastp', color_threshold=cut_height)
    plt.savefig(os.path.join(dirpath, rmsd_mat_obj_fname + '_clust_' + cluster_mode + '_dendrogram.png'))
    plt.close()



def manual_cluster(rmsd_mat_obj, rmsd_sq_mat, linkage, dirpath, num_clusters=2, cluster_mode='manual', cut_height=None):

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
            rmsd_val = rmsd_sq_mat[ind_pair[0], ind_pair[1]]
            rmsd_clust[num].append(rmsd_val)

    generate_cluster_data(clusters, rmsd_clust, rmsd_mat_obj, dirpath, cluster_mode=cluster_mode)
    cluster_rmsd_avg = np.mean(np.concatenate(rmsd_clust))
    print('Cluster rmsd average: ', cluster_rmsd_avg)
    dist_range, rmsd_avg_all, rmsd_avg_between_clust = dist_to_rmsd(linkage, rmsd_mat_sqform, rmsd_dist_mat, dirpath,
                                                                    cluster_mode=cluster_mode)

    if cut_height == None:
        cut_height = cut_ht_for_optimal_clusters(dist_range, rmsd_avg_all, cluster_rmsd_avg, upper_lim=-7)
        print('Cut height: ', cut_height)
    else:
        cut_height = cut_height

    # plot figure

    plt.figure()
    dn = hclust.dendrogram(linkage, distance_sort='ascending', truncate_mode='lastp', color_threshold=cut_height)
    plt.savefig(os.path.join(dirpath, rmsd_mat_obj_fname + '_clust_' + cluster_mode + '_dendrogram.png'))
    plt.close()



def f_test_cluster_numbers(linkage, rmsd_sq_mat, dirpath, max_num_clusters=10):
    min_num_clusters=4
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



if __name__ == '__main__':
    rmsd_mat_obj_fname = 'rmsd_mat_1.obj'
    dirpath = r"T:\Sugyan\Chae_RMSD"
    outdir = r"T:\Sugyan\Chae_RMSD\_1"
    rmsd_mat_obj_fpath = os.path.join(dirpath, rmsd_mat_obj_fname)

    rmsd_mat_obj = load_pickle_object(rmsd_mat_obj_fpath)
    rmsd_mat = rmsd_mat_obj.rmsd_list
    rmsd_mat_sqform = squareform(rmsd_mat)

    plot_rmsd_dist_mat_scatter(rmsd_mat_sqform, outdir)

    plot_rmsd_dist_hist(rmsd_mat, outdir)

    rmsd_dist_mat = pdist(rmsd_mat_sqform, 'euclidean')

    linkage = hclust.linkage(rmsd_dist_mat, method='average')

    autogen_cluster(rmsd_mat_obj, rmsd_mat_sqform, linkage, outdir, max_num_clusters=10)
    # manual_cluster(rmsd_mat_obj, rmsd_mat_sqform, linkage, outdir, num_clusters=5, cut_height=7)

    # num_cluster, clusters, cluster_rmsd = f_test_cluster_numbers(linkage, rmsd_mat_sqform, outdir, 10)
    # print ('Number of clusters: ', num_cluster)
    # generate_cluster_data(clusters, cluster_rmsd, rmsd_mat_obj, outdir)
    #
    # cluster_rmsd_avg = np.mean(np.concatenate(cluster_rmsd))
    # print('Cluster rmsd average: ', cluster_rmsd_avg)
    # dist_range, rmsd_avg_all, rmsd_avg_between_clust = dist_to_rmsd(linkage, rmsd_mat_sqform, rmsd_dist_mat, outdir)
    # cut_height = cut_ht_for_optimal_clusters(dist_range, rmsd_avg_all, cluster_rmsd_avg)
    # print('Cut height: ', cut_height)

    # cut_height = 200
    # plt.figure()
    # dn = hclust.dendrogram(linkage, distance_sort='ascending', truncate_mode='lastp', color_threshold=cut_height)
    # plt.savefig(os.path.join(outdir, rmsd_mat_obj_fname+'_clust_dendrogram.png'))
    # plt.close()

    print('heho')





    # c, coph_dists = cophenet(linkage, rmsd_dist_mat)

