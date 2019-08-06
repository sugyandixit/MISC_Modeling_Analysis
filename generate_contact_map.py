import os
import numpy as np
import matplotlib.pyplot as plt

def contact_map(dist_mat, dist_cut_off):
    contact_map_arr = np.zeros(np.shape(dist_mat))
    for ind, arr in enumerate(dist_mat):
        for ind_j, val in enumerate(arr):
            if type(dist_cut_off) == int:
                if val > dist_cut_off:
                    contact_map_arr[ind, ind_j] = 0
                elif val == 0:
                    contact_map_arr[ind, ind_j] = 0
                else:
                    contact_map_arr[ind, ind_j] = 1
            else:
                if val > dist_cut_off[1]:
                    contact_map_arr[ind, ind_j] = 0
                elif val < dist_cut_off[0]:
                    contact_map_arr[ind, ind_j] = 0
                elif val == 0:
                    contact_map_arr[ind, ind_j] = 0
                else:
                    contact_map_arr[ind, ind_j] = 1
    return contact_map_arr


def plot_contact_map(contact_map, dist_cut_off, pdbname, dirpath):
    plt.imshow(contact_map, cmap='binary', interpolation='nearest', origin='lower')
    # plt.show()
    plt.savefig(os.path.join(dirpath, pdbname+'_ct_map_' + str(dist_cut_off)+ '_.png'))
    plt.close()