# group contact maps based on cluster

import os
import numpy as np
import pandas as pd
import pickle
import generate_contact_map_list_object
from generate_contact_map_list_object import ContactMap


class ClusterContactMaps(object):
    """
    store contact maps by cluster
    """
    def __init__(self):
        self.cluster_pdb_index_num = []
        self.cluster_contact_maps = []



def save_object_to_pickle(obj, file_path):
    with open(file_path, 'wb') as out_obj:
        pickle.dump(obj, out_obj)
        out_obj.close()


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

def generate_cluster_contact_maps(dataframe, contact_mat_object):
    cluster_contact_map_obj = ClusterContactMaps()
    cluster_id = np.unique(dataframe['clust_num'].values)
    cluster_contact_arr = [[] for i in range(len(cluster_id))]
    cluster_pdb_ind_arr = [[] for i in range(len(cluster_id))]
    for num, clust_id in enumerate(cluster_id):
        clust_df = dataframe[dataframe['clust_num']==clust_id]
        for index, clust_pdb_fpath in enumerate(clust_df['pdb_fname']):
            clust_pdb_name = os.path.split(clust_pdb_fpath)[1]
            for ind, (contact_mat_pdb_name, contact_mat_) in enumerate(zip(contact_mat_object.pdbname, contact_mat_object.contact_mat)):
                if contact_mat_pdb_name == clust_pdb_name:
                    cluster_contact_arr[num].append(contact_mat_)
                    cluster_pdb_ind_arr[num].append(ind)
                else:
                    pass
    cluster_contact_map_obj.cluster_contact_maps = cluster_contact_arr
    cluster_contact_map_obj.cluster_pdb_index_num = cluster_pdb_ind_arr
    return cluster_contact_map_obj


if __name__=='__main__':
    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\cluster_1\num_clusters_6"
    contact_mat_obj_fname = "contact_mat_5.obj"
    cluster_pdb_fname = "cluster_pdb_files.csv"
    df = pd.read_csv(os.path.join(dirpath, cluster_pdb_fname))
    df = df.drop_duplicates(subset='pdb_fname')
    contact_mat_obj = load_pickle_object(os.path.join(dirpath, contact_mat_obj_fname))
    cluster_ct_map_obj = generate_cluster_contact_maps(df, contact_mat_obj)
    cluster_ct_map_obj_fname = 'cluster_'+contact_mat_obj_fname
    save_object_to_pickle(cluster_ct_map_obj, os.path.join(dirpath, cluster_ct_map_obj_fname))