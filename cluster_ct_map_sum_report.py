# sum contact maps in each cluster to get overall frequency matrix and plot imshow and bar graph

import os
import pickle
import numpy as np
import pandas as pd
from scipy.stats import zscore
import matplotlib.pyplot as plt
import cluster_contact_maps
from cluster_contact_maps import ClusterContactMaps
from generate_contact_map_list_object import ContactMap


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

def plot_imshow(array, seq_xaxis, seq_yaxis, fname, dirpath):
    # plt.figure(figsize=(8,10))
    plt.imshow(array, origin='lower')
    plt.xticks(np.arange(0, len([x for x in seq_xaxis])), seq_xaxis)
    plt.yticks(np.arange(0, len([x for x in seq_yaxis]), 2), seq_yaxis)
    plt.colorbar()
    plt.savefig(os.path.join(dirpath, fname))
    plt.close()

def plot_bar_graph(array, seq, fname, dirpath, color='red'):
    ind = np.arange(1, len(array)+1)
    plt.bar(ind, array, width=1, color=color)
    plt.xticks(ind, seq)
    plt.savefig(os.path.join(dirpath, fname))
    plt.close()


def write_zscore_output(resseq, zscore_arr, fname, dirpath):
    outlines = ''
    header = '{},{},{}\n'.format('residue_num', 'residue', 'zscore')
    outlines += header
    for ind, (resid, zscore) in enumerate(zip(resseq, zscore_arr)):
        line = '{},{},{}\n'.format(ind+1, resid, zscore)
        outlines += line

    with open(os.path.join(dirpath, fname), 'w') as zscore_out:
        zscore_out.write(outlines)
        zscore_out.close()

    print('test')


def write_zscore_for_pairs(z_score_pairs, ind_1, ind_2, res_sq_1, res_sq_2, outid, dirpath):

    output_string = ''
    header = 'res1_num,res2_num,res1_,res2_,zscore\n'
    output_string += header

    for ind, (res1num, res2num, zscore_) in enumerate(zip(ind_1, ind_2, z_score_pairs)):

        res1_ = res_sq_1[res1num]
        res2_ = res_sq_2[res2num]

        line = '{},{},{},{},{}\n'.format(res1num, res2num, res1_, res2_, zscore_)

        output_string += line


    out_name = 'pair_zscore_'+ outid + '.csv'
    with open(os.path.join(dirpath, out_name), 'w') as outfile:
        outfile.write(output_string)
        outfile.close()

def construct_grid(wh, wv, data):
    wv_unique = np.unique(wv)
    wh_unique = np.unique(wh)
    data_grid = np.zeros(shape=(len(wv_unique), len(wh_unique)))
    ind_new = 0
    for ind, (wv_, wh_) in enumerate(zip(wv, wh)):
        for i, wv_un in enumerate(wv_unique):
            for j, wh_un in enumerate(wh_unique):
                if wv_un == wv_:
                    if wh_un == wh_:
                        data_grid[i, j] = data[ind_new]
                        ind_new += 1

    data_grid[data_grid == 0] = 'nan'
    return data_grid

def plot_zscore_imhow(zscoremat, cmap, distcutoff, xaxis, yaxis, dirpath):
    y_ticks = np.arange(0, len(yaxis), 2)
    y_tick_labels = []
    for ind in y_ticks:
        y_tick_labels.append(yaxis[ind])
    plt.imshow(zscoremat, origin=('lower', 'lower'), aspect='auto', cmap=cmap)
    plt.xticks(np.arange(0, len(xaxis), 1), xaxis)
    plt.yticks(y_ticks, y_tick_labels)
    plt.colorbar()
    plt.savefig(os.path.join(dirpath, 'pair_zscore_imshow_dist_cut_off_'+str(distcutoff)+'.pdf'))
    plt.close()


def z_score_for_pairs(sum_ct, resseq1, resseq2, dirpath, outid, image_type):

    ind_1 = []
    ind_2 = []
    sum_ct_arr = []
    for ind1 in range(len(sum_ct[0,:])):
        for ind2 in range(len(sum_ct[:, 0])):
            ind_1.append(ind1)
            ind_2.append(ind2)
            sum_ct_val = sum_ct[ind2, ind1]
            sum_ct_arr.append(sum_ct_val)

    sum_ct_arr = np.array(sum_ct_arr)
    z_score_pairs = zscore(sum_ct_arr)
    sort_ind = np.argsort(z_score_pairs)
    z_score_pair_sort = z_score_pairs[sort_ind[::-1]]
    ind_1_sort = np.array(ind_1)[sort_ind[::-1]]
    ind_2_sort = np.array(ind_2)[sort_ind[::-1]]

    write_zscore_for_pairs(z_score_pair_sort, ind_1_sort, ind_2_sort, resseq1, resseq2, outid, dirpath)

    z_score_sq_mat = construct_grid(ind_1_sort, ind_2_sort, z_score_pair_sort)

    z_score_mat_masked = np.ma.masked_where(z_score_sq_mat <= 2, z_score_sq_mat)
    # mask_array = np.ma.array(z_score_mat, mask=np.where(z_score_mat <=2))
    cmap = plt.cm.get_cmap('PuBu')
    cmap.set_bad(color='white')

    plot_zscore_imhow(z_score_mat_masked, cmap, distcutoff=outid, xaxis=resseq1, yaxis=resseq2, dirpath=dirpath)

    # z_score_sq_mat = z_score_pairs.reshape(len(sum_ct[0,:]), len(sum_ct[:,0]))

    # plot_imshow(z_score_sq_mat, resseq1, resseq2, 'pair_zscore_imshow_'+outid + image_type, dirpath)

    return z_score_sq_mat





    print('heho')




def total_contacts_for_cluster(contact_map_obj, resseq1, resseq2, dirpath, outid, image_type='.png'):
    """

    :param contact_map_obj:
    :param resseq1:
    :param resseq2:
    :param dirpath:
    :param image_type:
    :return:
    """

    sum_ct = np.sum(contact_map_obj.contact_mat, axis=0)

    z_score_for_pairs(sum_ct, resseq1, resseq2, dirpath, outid, image_type)

    plot_imshow(sum_ct, resseq1, resseq2, 'total_contacts_cluster_'+outid+image_type, dirpath)

    ct_dataframe = pd.DataFrame(sum_ct)
    ct_dataframe.to_csv(os.path.join(dirpath, 'total_contacts_cluster_' + outid + '.csv'))

    mol1_freq = np.sum(sum_ct, axis=0)
    mol2_freq = np.sum(sum_ct, axis=1)

    # plot_bar_graph(mol1_freq, resseq1, 'mol1_total_contacts_cluster_' + outid + image_type, dirpath, color='blue')
    # plot_bar_graph(mol2_freq, resseq2, 'mol2_total_contacts_cluster_' + outid + image_type, dirpath, color='red')

    mol1_zscore = zscore(mol1_freq)
    mol2_zscore = zscore(mol2_freq)

    write_zscore_output(resseq1, mol1_zscore, 'mol1_zscore_cluster_' + outid + '.csv', dirpath)
    write_zscore_output(resseq2, mol2_zscore, 'mol2_zscore_cluster_' + outid + '.csv', dirpath)

    # plot_bar_graph(mol1_zscore, resseq1, 'mol1_total_contacts_zscore_cluster_' + outid + image_type, dirpath,
    #                color='blue')
    # plot_bar_graph(mol2_zscore, resseq2, 'mol2_total_contacts_zscore_cluster_' + outid + image_type, dirpath,
    #                color='red')




def total_contacts_for_cluster_using_combined_cluster_ct_object(cluster_ct_object, resseq1, resseq2, dirpath, image_type='.png'):


    for index, cluster_ct in enumerate(cluster_ct_object.cluster_contact_maps):

        sum_ct = np.sum(cluster_ct, axis=0)

        plot_imshow(sum_ct, resseq1, resseq2, 'total_contacts_cluster_'+str(index)+image_type, dirpath)

        z_score_for_pairs(sum_ct, resseq1, resseq2, dirpath, outid='clust_'+str(index), image_type=image_type)

        ct_dataframe = pd.DataFrame(sum_ct)
        ct_dataframe.to_csv(os.path.join(dirpath, 'total_contacts_cluster_' + str(index) + '.csv'))

        mol1_freq = np.sum(sum_ct, axis=0)
        mol2_freq = np.sum(sum_ct, axis=1)

        plot_bar_graph(mol1_freq, resseq1, 'mol1_total_contacts_cluster_'+str(index)+image_type, dirpath, color='blue')
        plot_bar_graph(mol2_freq, resseq2, 'mol2_total_contacts_cluster_' + str(index) + image_type, dirpath, color='red')

        mol1_zscore = zscore(mol1_freq)
        mol2_zscore = zscore(mol2_freq)

        write_zscore_output(resseq1, mol1_zscore, 'mol1_zscore_cluster_'+str(index)+'.csv', dirpath)
        write_zscore_output(resseq2, mol2_zscore, 'mol2_zscore_cluster_' + str(index) + '.csv', dirpath)

        plot_bar_graph(mol1_zscore, resseq1, 'mol1_total_contacts_zscore_cluster_' + str(index) + image_type, dirpath, color='blue')
        plot_bar_graph(mol2_zscore, resseq2, 'mol2_total_contacts_zscore_cluster_' + str(index) + image_type, dirpath, color='red')



if __name__=='__main__':
    ## for serf ab interaction analysis
    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\cluster_5"
    cluster_ct_obj_fname = "cluster_contact_mat_5.obj"
    mol1_res = 'DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVV'
    mol2_res = 'MARGNQRDLARQKNLKKQKDMAKNQKKSGDPKKRMESDAEILRQKQAAADARREAEKLEK'
    cluster_ct_obj = load_pickle_object(os.path.join(dirpath, cluster_ct_obj_fname))
    total_contacts_for_cluster_using_combined_cluster_ct_object(cluster_ct_obj, mol1_res, mol2_res, dirpath)

    #for single cluster object analysis
    # dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\rmsd_mat\cluster_manual_7"
    # cluster_ct_obj_fname = "contact_mat_clust_7_[5, 10].obj"
    # mol1_res = 'MARGNQRDLARQKNLKKQKDMAKNQKKSGDPKKRMESDAEILRQKQAAADARREAEKLEKLKAEKTRR'
    # mol2_res = mol1_res
    # clust_ct_obj = load_pickle_object(os.path.join(dirpath, cluster_ct_obj_fname))
    # total_contacts_for_cluster(clust_ct_obj, mol1_res, mol2_res, dirpath, outid='cluster_man_7_[5, 10]', image_type='.pdf')


