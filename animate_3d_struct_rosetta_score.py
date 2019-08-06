# animate the landscape


import os
import pandas as pd
import pickle
import plot_landscape_3d
from plot_landscape_3d import plot_landscape, generate_grid_data
import numpy as np
import rmsd_mat_obj
from rmsd_mat_obj import RMSDMAT
from scipy.spatial.distance import squareform, pdist
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation


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

    rmsd_mat_obj_ = load_pickle_object(os.path.join(dirpath, rmsd_mat_obj_fname))

    rmsd_mat = rmsd_mat_obj_.rmsd_list
    rmsd_mat_sqform = squareform(rmsd_mat)

    rmsd_dist_mat = pdist(rmsd_mat_sqform, 'euclidean')

    similarity = squareform(rmsd_dist_mat)

    grid_data = generate_grid_data(similarity, ros_score, dirpath)

    fig = plt.figure()
    ax = Axes3D(fig)

    def init():
        ax.plot_surface(grid_data[0], grid_data[1], grid_data[2], rstride=1, cstride=1, cmap='jet')
        return ()

    plt.axis('off')

    def animate(i):
        ax.view_init(elev=13, azim=i+60)
        return()

    anim = animation.FuncAnimation(fig, animate, init_func=init, frames=720, interval=10, blit=True)
    anim.save(os.path.join(dirpath, 'test_animation.mp4'), fps=30)