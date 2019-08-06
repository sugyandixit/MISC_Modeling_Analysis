import os
import pickle
import numpy as np
# import distance_map
# from distance_map import generate_distance_map
from generate_distance_map_list_object import DistanceMap
import generate_contact_map
from generate_contact_map import contact_map, plot_contact_map
# from dataclasses import dataclass

# @dataclass
class ContactMap(object):

    def __init__(self):

        self.dirpath = None
        self.pdbname = None
        self.chain_id = None
        self.dist_cut_off = None
        self.contact_mat = None


def save_object_to_pickle(obj, file_path):
    with open(file_path, 'wb') as out_obj:
        pickle.dump(obj, out_obj)
        out_obj.close()


def get_distance_matrix_obj_files(dirpath):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('distance_matrix.obj')]
    return files


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

def generate_contact_mat_list_object(dist_mat_obj_fpath, dist_cut_off, out_id, plot_contacts=False):
    """
    generate contact mat list object and save
    :param dirpath: directory to save cotnact mat list obj
    :param dist_mat_obj_fpath: distance matrix object fpath
    :param dist_cut_off: distance cutoff used to generate contacts
    :param out_id: out id for saving file
    :return: contact mat list object
    """
    dist_mat_obj_dirpath = os.path.split(dist_mat_obj_fpath)[0]
    dist_mat = load_pickle_object(dist_mat_obj_fpath)
    contact_map_list = []
    dir_list = []
    pdb_list = []
    chain_id_list = []
    for ind, (dist_matrix, dirpath, pdbname, chain_id) in enumerate(
            zip(dist_mat.distance_mat, dist_mat.dirpath, dist_mat.pdbname, dist_mat.chain_id)):
        ct_map = contact_map(dist_matrix, dist_cut_off)
        if np.any(ct_map > 0):
            print(str(dist_cut_off) + 'Angstrom Contact matrix for pdb --> ' + pdbname)
            contact_map_list.append(ct_map)
            dir_list.append(dirpath)
            pdb_list.append(pdbname)
            chain_id_list.append(chain_id)
            if plot_contacts:
                plot_contact_map(ct_map, dist_cut_off, pdbname, dirpath)

    ct_map_obj = ContactMap()
    ct_map_obj.dist_cut_off = dist_cut_off
    ct_map_obj.dirpath = dir_list
    ct_map_obj.pdbname = pdb_list
    ct_map_obj.chain_id = chain_id_list
    ct_map_obj.contact_mat = contact_map_list

    ct_map_obj_fpath = os.path.join(dist_mat_obj_dirpath, 'contact_mat_' + out_id + '_' + str(dist_cut_off) + '.obj')
    save_object_to_pickle(ct_map_obj, ct_map_obj_fpath)

    return ct_map_obj



if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\rmsd_mat\cluster_manual_7"
    distance_mat_files = get_distance_matrix_obj_files(dirpath)

    for ind, file in enumerate(distance_mat_files):

        distance_mat_fpath = os.path.join(dirpath, file)
        out_id_ = 'clust_7'

        generate_contact_mat_list_object(distance_mat_fpath, dist_cut_off=[5,10], out_id=out_id_, plot_contacts=False)