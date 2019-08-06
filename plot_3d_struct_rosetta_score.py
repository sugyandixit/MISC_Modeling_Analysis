# plot 3d landscape with rosetta energy score

import os
import pandas as pd
import pickle
import plot_landscape_3d
from plot_landscape_3d import plot_landscape
import numpy as np
import rmsd_mat_obj
from rmsd_mat_obj import RMSDMAT
from scipy.spatial.distance import squareform, pdist


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

if __name__=='__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models"
    rmsd_mat_obj_fname = 'rmsd_mat_ct_dist_5.obj'
    pdb_ros_score_fname =  "pdb_list_ct_map_5.csv"

    df = pd.read_csv(os.path.join(dirpath, pdb_ros_score_fname), sep=',')
    ros_score = df['rosetta_score'].values

    rmsd_mat_obj = load_pickle_object(os.path.join(dirpath, rmsd_mat_obj_fname))

    rmsd_mat = rmsd_mat_obj.rmsd_list
    rmsd_mat_sqform = squareform(rmsd_mat)

    rmsd_dist_mat = pdist(rmsd_mat_sqform, 'euclidean')

    similarity = squareform(rmsd_dist_mat)

    plot_landscape(similarity, ros_score, dirpath)

