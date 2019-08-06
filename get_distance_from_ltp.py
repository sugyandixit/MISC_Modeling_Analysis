import os
import numpy as np
import pickle
import pdb_utility
from pdb_utility import ParsePDB, distance_between_two_points
import matplotlib.pyplot as plt
from cluster_ct_map_sum_report import z_score_for_pairs
from generate_contact_map import contact_map


def get_files(dirpath, startid, endid):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith(endid)]
    if startid:
        files = [x for x in files if x.startswith(startid)]
    return files


def get_dist_back_forth(pdb_read_obj, ref_res='LTP', ref_res_atom_name='NZ', res_atom_name='O', res_span=[-5,5]):
    """
    get distance from
    :param pdb_read_obj: pdb read object
    :param ref_res: reference residue name
    :param ref_res_atom_name: reference residue atom name
    :param res_atom_name: restidue atom name . distance is calcualted between ref rest atom name and rest atom name
    :param res_span: residues to span from ref residue
    :return: distance list
    """

    res_span_ = np.arange(res_span[0], res_span[1]+1, 1)

    ref_res_num_arr = []
    ref_res_name_arr = []
    ref_coor_arr = []

    for index1, (pdb_res_name, pdb_atom_name) in enumerate(zip(pdb_read_obj.atm_resname, pdb_read_obj.atm_name)):
        if pdb_res_name == ref_res:
            if pdb_atom_name == ref_res_atom_name:
                ref_res_num = int(pdb_read_obj.atm_resseqnum[index1])
                ref_res_name = pdb_read_obj.atm_resname[index1]
                ref_atom_coords = [pdb_read_obj.atm_xcoor[index1],
                                   pdb_read_obj.atm_ycoor[index1],
                                   pdb_read_obj.atm_zcoor[index1]]

                ref_res_num_arr.append(ref_res_num)
                ref_res_name_arr.append(ref_res_name)
                ref_coor_arr.append(ref_atom_coords)


    pdb_res_num_array = np.array(pdb_read_obj.atm_resseqnum, dtype='int')
    min_res_num = pdb_res_num_array.min()
    max_res_num = pdb_res_num_array.max()

    dist_list_array = np.zeros((len(ref_res_num_arr), len(res_span_)))


    for index2, (ref_res_num_, ref_coor_) in enumerate(zip(ref_res_num_arr, ref_coor_arr)):
        res_span_list = ref_res_num_ + res_span_
        for index3, res_span_atom_num in enumerate(res_span_list):
            if res_span_atom_num >= min_res_num:
                if res_span_atom_num <= max_res_num:
                    for index4, (pdb_res_num, pdb_atom_name) in enumerate(zip(pdb_read_obj.atm_resseqnum, pdb_read_obj.atm_name)):
                        if int(pdb_res_num) == res_span_atom_num:
                            if pdb_atom_name == res_atom_name:
                                res_coor = [pdb_read_obj.atm_xcoor[index4],
                                            pdb_read_obj.atm_ycoor[index4],
                                            pdb_read_obj.atm_zcoor[index4]]
                                dist_ = distance_between_two_points(ref_coor_, res_coor)
                                dist_list_array[index2][index3] = dist_


    dist_dict = dict()
    dist_dict['ref_res'] = ref_res
    dist_dict['ref_res_atom_name'] = ref_res_atom_name
    dist_dict['res_atom_name'] = res_atom_name
    dist_dict['res_span'] = res_span_
    dist_dict['ref_res_num_arr'] = ref_res_num_arr
    dist_dict['ref_res_name_arr'] = ref_res_name_arr
    dist_dict['dist_list_array'] = dist_list_array

    return dist_dict


def get_dist_one_to_all(pdb_read_obj, ref_res='LTP', ref_res_atom_name='NZ', res_atom_name='O'):
    """
    get distance from ref res atom name to all other res atom name
    :param pdb_read_obj: pdb read object
    :param ref_res: ref res
    :param ref_res_atom_name: ref rest atom name
    :param res_atom_name: rest aof res atom name
    :return: dict
    """

    ref_res_num_arr = []
    ref_res_name_arr = []
    ref_coor_arr = []

    for index1, (pdb_res_name, pdb_atom_name) in enumerate(zip(pdb_read_obj.atm_resname, pdb_read_obj.atm_name)):
        if pdb_res_name == ref_res:
            if pdb_atom_name == ref_res_atom_name:
                ref_res_num = int(pdb_read_obj.atm_resseqnum[index1])
                ref_res_name = pdb_read_obj.atm_resname[index1]
                ref_atom_coords = [pdb_read_obj.atm_xcoor[index1],
                                   pdb_read_obj.atm_ycoor[index1],
                                   pdb_read_obj.atm_zcoor[index1]]

                ref_res_num_arr.append(ref_res_num)
                ref_res_name_arr.append(ref_res_name)
                ref_coor_arr.append(ref_atom_coords)

    pdb_res_num_array = np.array(pdb_read_obj.atm_resseqnum, dtype='int')
    uniq_res_num_arr = np.unique(pdb_res_num_array)
    dist_list_array = np.zeros((len(ref_res_num_arr), len(uniq_res_num_arr)))



    for index2, (ref_res_num_, ref_coor_) in enumerate(zip(ref_res_num_arr, ref_coor_arr)):
        for index4, (pdb_res_num, pdb_atom_name) in enumerate(zip(pdb_read_obj.atm_resseqnum, pdb_read_obj.atm_name)):
            if pdb_atom_name == res_atom_name:
                res_num = int(pdb_res_num)
                res_coor = [pdb_read_obj.atm_xcoor[index4],
                            pdb_read_obj.atm_ycoor[index4],
                            pdb_read_obj.atm_zcoor[index4]]
                dist_ = distance_between_two_points(ref_coor_, res_coor)
                dist_list_array[index2][res_num-1] = dist_


    dist_dict = dict()
    dist_dict['ref_res'] = ref_res
    dist_dict['ref_res_atom_name'] = ref_res_atom_name
    dist_dict['res_atom_name'] = res_atom_name
    dist_dict['ref_res_num_arr'] = ref_res_num_arr
    dist_dict['ref_res_name_arr'] = ref_res_name_arr
    dist_dict['res_num_arr'] = uniq_res_num_arr
    dist_dict['dist_list_array'] = dist_list_array

    return dist_dict

def plot_dist_contour(dist_mat, dirpath, outid, xaxis, yaxis):
    y_ticks = np.arange(0, len(yaxis), 2)
    y_tick_labels = []
    for ind in y_ticks:
        y_tick_labels.append(yaxis[ind])
    plt.imshow(dist_mat, origin=('lower', 'lower'), aspect='auto')
    plt.xticks(np.arange(0, len(xaxis), 1), xaxis)
    plt.yticks(y_ticks, y_tick_labels)
    plt.colorbar()
    plt.savefig(os.path.join(dirpath, outid + '_dist_map.pdf'))
    plt.close()

# def contact_map(dist_mat, dist_cut_off):
#     contact_map_arr = np.zeros(np.shape(dist_mat))
#     for ind, arr in enumerate(dist_mat):
#         for ind_j, val in enumerate(arr):
#             if val > dist_cut_off:
#                 contact_map_arr[ind, ind_j] = 0
#             elif val == 0:
#                 contact_map_arr[ind, ind_j] = 0
#             else:
#                 contact_map_arr[ind, ind_j] = 1
#     return contact_map_arr

def plot_contact_map(contact_map, dist_cut_off, pdbname, dirpath, xaxis, yaxis):
    y_ticks = np.arange(0, len(yaxis), 2)
    y_tick_labels = []
    for ind in y_ticks:
        y_tick_labels.append(yaxis[ind])
    plt.imshow(contact_map, cmap='binary', interpolation='nearest', origin='lower', aspect='auto')
    # plt.show()
    plt.xticks(np.arange(0, len([x for x in xaxis])), xaxis)
    plt.yticks(y_ticks, y_tick_labels)
    plt.savefig(os.path.join(dirpath, pdbname+'_ct_map_' + str(dist_cut_off)+ '_.pdf'))
    plt.close()

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


def write_dist_output_one_to_all_from_list_of_pdb_files(dirpath, list_pdb_files, outid='dist_one_to_all_list_output', ref_res='LTP', ref_res_atom_name='NZ', res_atom_name='O', dist_cut_off=4):



    dist_arr_tp_arr = []

    ref_res_num_arr = []
    res_num_arr = []

    for index, file in enumerate(list_pdb_files):
        print('Calculating distances for pdb_file ', str(file))
        pdb_obj = ParsePDB(file, dirpath).read()
        dist_dict = get_dist_one_to_all(pdb_obj, ref_res=ref_res, ref_res_atom_name=ref_res_atom_name, res_atom_name=res_atom_name)

        ref_res_num_arr = dist_dict['ref_res_num_arr']
        res_num_arr = dist_dict['res_num_arr']

        dist_arr_tp = dist_dict['dist_list_array'].T
        dist_arr_tp_arr.append(dist_arr_tp)


    contact_arr = get_contact_arr(dist_arr_tp_arr, dist_cut_off)

    sum_contact_arr = np.sum(contact_arr, axis=0)

    z_score_mat = z_score_for_pairs(sum_contact_arr, ref_res_num_arr, res_num_arr, dirpath, 'dist_cut_off_'+str(dist_cut_off), '.pdf')

    z_score_mat_masked = np.ma.masked_where(z_score_mat <= 1, z_score_mat)
    # mask_array = np.ma.array(z_score_mat, mask=np.where(z_score_mat <=2))
    cmap = plt.cm.get_cmap('PuBu')
    cmap.set_bad(color='white')

    plot_zscore_imhow(z_score_mat_masked, cmap, dist_cut_off, ref_res_num_arr, res_num_arr, dirpath)



    mean_dist_arr = np.mean(dist_arr_tp_arr, axis=0)
    std_dist_arr = np.std(dist_arr_tp_arr, axis=0)

    dist_header = ','.join(str(x) for x in ref_res_num_arr)
    header = 'res_num,'+dist_header+'\n'

    data_string = ''
    for index, (dist_) in enumerate(mean_dist_arr):
        dist_str = ','.join([str(x) for x in dist_])
        line = '{},{}\n'.format(index+1,
                                   dist_str)
        data_string += line
    output_string = header + data_string
    outname = outid + '_mean.csv'
    with open(os.path.join(dirpath, outname), 'w') as outfile:
        outfile.write(output_string)
        outfile.close()

    data_string = ''
    for index, (dist_) in enumerate(std_dist_arr):
        dist_str = ','.join([str(x) for x in dist_])
        line = '{},{}\n'.format(index + 1,
                                dist_str)
        data_string += line
    output_string = header + data_string
    outname = outid + '_std.csv'
    with open(os.path.join(dirpath, outname), 'w') as outfile:
        outfile.write(output_string)
        outfile.close()



    plot_dist_contour(mean_dist_arr, dirpath, outid='mean', xaxis=ref_res_num_arr, yaxis=res_num_arr)
    plot_dist_contour(std_dist_arr, dirpath, outid='std', xaxis=ref_res_num_arr, yaxis=res_num_arr)

    plot_contact_map(contact_map(mean_dist_arr, dist_cut_off), dist_cut_off, 'mean', dirpath, xaxis=ref_res_num_arr,
                     yaxis=res_num_arr)



def get_contact_arr(list_of_dist_arr, dist_cut_off):

    contact_arr = []

    for dist_arr in list_of_dist_arr:

        contact_map_arr = contact_map(dist_arr, dist_cut_off)

        contact_arr.append(contact_map_arr)

    return contact_arr








def write_dist_output_back_forth_from_list_of_pdb_files(dirpath, list_pdb_files, outid='dist_list_output', ref_res='LTP', ref_res_atom_name='NZ', res_atom_name='O', res_span=[-5,5]):

    res_span_atoms = np.arange(res_span[0], res_span[1] + 1, 1)
    # res_span_atom_str = ','.join([str(x) for x in res_span_atoms])
    res_span_atom_str = []
    for res_span_num in res_span_atoms:
        if res_span_num < 0:
            span_str = 'i-'+str(abs(res_span_num))
        if res_span_num > 0:
            span_str = 'i+'+str(res_span_num)
        if res_span_num == 0:
            span_str = 'i+0'
        res_span_atom_str.append(span_str)
    res_span_atom_str = ','.join([x for x in res_span_atom_str])

    header = 'pdb_file,ref_res,ref_res_num,ref_res_atom_name,res_atom_name,'+res_span_atom_str+'\n'

    data_string = ''

    for index, file in enumerate(list_pdb_files):

        print('Calculating distances for pdb_file ', str(file))

        pdb_obj = ParsePDB(file, dirpath).read()
        dist_dict = get_dist_back_forth(pdb_obj, ref_res=ref_res, ref_res_atom_name=ref_res_atom_name, res_atom_name=res_atom_name, res_span=res_span)

        for index2, (ref_res_num, ref_res_name, dist_list) in enumerate(zip(dist_dict['ref_res_num_arr'], dist_dict['ref_res_name_arr'], dist_dict['dist_list_array'])):

            dist_str = ','.join([str(x) for x in dist_list])

            line = '{},{},{},{},{},{}\n'.format(file,
                                                ref_res_name,
                                                ref_res_num,
                                                ref_res_atom_name,
                                                res_atom_name,
                                                dist_str)

            data_string += line

    output_string = header + data_string

    outname = outid + '.csv'

    with open(os.path.join(dirpath, outname), 'w') as outfile:

        outfile.write(output_string)
        outfile.close()



if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_nomod\rmsdmat\cluster_manual_7"
    # pdb_name = 'repid_0_frame_12101.pdb'
    # pdb_obj = ParsePDB(pdb_name, dirpath).read()

    # get_dist_one_to_all(pdb_obj)

    pdb_file_list = get_files(dirpath, startid='repid', endid='.pdb')
    write_dist_output_one_to_all_from_list_of_pdb_files(dirpath, pdb_file_list, ref_res='LYS', dist_cut_off=[5,10])

    # write_dist_output_back_forth_from_list_of_pdb_files(dirpath, pdb_file_list)