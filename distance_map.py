import os
import itertools
import numpy as np
import matplotlib.pyplot as plt
from pdb_utility import ParsePDB, distance_between_two_points

def combination_array(list1, list2):
    comb_arr = []
    for item1 in list1:
        for item2 in list2:
            comb_arr.append([item1, item2])
    return comb_arr


def distance_list(coor_comb_arr):
    dist_list = []
    for index in range(len(coor_comb_arr)):
        coor1 = coor_comb_arr[index][0]
        coor2 = coor_comb_arr[index][1]
        distance = distance_between_two_points(coor1, coor2)
        dist_list.append(distance)
    return dist_list


def distance_mat(res_comb_arr, dist_list):
    id_1 = []
    id_2 = []
    for index, res_comb in enumerate(res_comb_arr):
        id_1.append(res_comb[0][2])
        id_2.append(res_comb[1][2])
    id_x = np.unique(id_1)
    id_y = np.unique(id_2).T
    xx, yy = np.meshgrid(id_x, id_y)
    # yy = yy.T
    distmat = np.zeros(np.shape(xx))

    for num1 in range(len(id_x)):
        for num2 in range(len(id_y)):
            for ind, (res_comb, dist) in enumerate(zip(res_comb_arr, dist_list)):
                if id_x[num1] == res_comb[0][2]:
                    if id_y[num2] == res_comb[1][2]:
                        distmat[num2, num1] = dist
    return distmat


def plot_contour_distmat(dist_mat, pdbname, dirpath):
    plt.imshow(dist_mat, origin=('lower','lower'), cmap='magma')
    plt.colorbar()
    plt.savefig(os.path.join(dirpath, pdbname + '_dist_map.png'))
    plt.close()


def write_dist_mat_out(res_comb_arr, dist_mat, dirpath):
    outline = ''
    header = 'chain_id_1,atm_name_1,res_seq_num_1,res_name_1,chain_id_2,atm_name_2,res_seq_num_2,res_name_2,dist\n'
    outline += header
    for index, (res_comb, dist) in enumerate(zip(res_comb_arr, dist_mat)):
        string_1 = ','.join(str(x) for x in res_comb[0])
        string_2 = ','.join(str(x) for x in res_comb[1])
        line = '{},{},{}\n'.format(string_1, string_2, str(dist))
        outline += line
    with open(os.path.join(dirpath, 'dist_mat_output.csv'), 'w') as outfile:
        outfile.write(outline)
        outfile.close()


def assemble_coordinate_combination_array(dirpath, pdbname, chain_ids=None):
    pdb_obj = ParsePDB(pdbname, dirpath).read()
    atom_name = 'CA'
    # chain_ids = ['A', 'S']
    if chain_ids:
        store_list = [[] for _ in range(len(chain_ids))]
        coors_list = [[] for _ in range(len(chain_ids))]

        for ind, chain_id in enumerate(chain_ids):
            for index in range(len(pdb_obj.atm_sernum)):
                if pdb_obj.atm_chainid[index] == chain_id:
                    if pdb_obj.atm_name[index] == atom_name:
                        chain_name = pdb_obj.atm_chainid[index]
                        atm_name = pdb_obj.atm_name[index]
                        res_name = pdb_obj.atm_resname[index]
                        res_num = int(pdb_obj.atm_resseqnum[index])
                        xcoor = float(pdb_obj.atm_xcoor[index])
                        ycoor = float(pdb_obj.atm_ycoor[index])
                        zcoor = float(pdb_obj.atm_zcoor[index])
                        store_list[ind].append([chain_name, atm_name, res_num, res_name])
                        coors_list[ind].append([xcoor, ycoor, zcoor])
    else:

        store_arr = []
        coors_arr = []

        for index in range(len(pdb_obj.atm_sernum)):
            if pdb_obj.atm_name[index] == atom_name:
                chain_name = 'J'
                atm_name = pdb_obj.atm_name[index]
                res_name = pdb_obj.atm_resname[index]
                res_num = int(pdb_obj.atm_resseqnum[index])
                xcoor = float(pdb_obj.atm_xcoor[index])
                ycoor = float(pdb_obj.atm_ycoor[index])
                zcoor = float(pdb_obj.atm_zcoor[index])
                store_arr.append([chain_name, atm_name, res_num, res_name])
                coors_arr.append([xcoor, ycoor, zcoor])

        store_list = [store_arr, store_arr]
        coors_list = [coors_arr, coors_arr]


    coor_comb_arr = combination_array(coors_list[0], coors_list[1])
    res_comb_arr = combination_array(store_list[0], store_list[1])
    return coor_comb_arr, res_comb_arr

def get_pdb_files(dirpath):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('.pdb')]
    return files


def generate_distance_map(pdbname, dirpath, chain_ids):
    coor_comb_arr, res_comb_arr = assemble_coordinate_combination_array(dirpath, pdbname, chain_ids)
    dist_list = distance_list(coor_comb_arr)
    dist_mat = distance_mat(res_comb_arr, dist_list)
    return dist_mat


if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\step2_remd_coor_final"
    pdb_files = get_pdb_files(dirpath)
    for index, pdbname in enumerate(pdb_files):
        dist_mat = generate_distance_map(pdbname, dirpath, chain_ids=None)
        plot_contour_distmat(dist_mat, pdbname, dirpath)
        print('heho')