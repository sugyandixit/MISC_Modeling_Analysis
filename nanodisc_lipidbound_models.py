import os
from pdb_utility import ParsePDB, ParseCRD, curate_pdb_segid, curate_crd_segid, convert_pdb_to_crd, write_pdb, write_crd, centerofmass, distance_between_two_points
import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

def compute_com_segid(readpdb_object, segid_name_list):
    xcoor = []
    ycoor = []
    zcoor = []
    for seg in segid_name_list:
        for ind, segid in enumerate(readpdb_object.atm_segid):
            if segid == seg:
                xcoor.append(readpdb_object.atm_xcoor[ind])
                ycoor.append(readpdb_object.atm_ycoor[ind])
                zcoor.append(readpdb_object.atm_zcoor[ind])
    com_segid = centerofmass(xcoor, ycoor, zcoor)
    return com_segid

def compute_distance_list(readpdb_obj, segid1_com, segid2_name):
    distance_list = []
    resnum_list = []
    for ind, segid in enumerate(readpdb_obj.atm_segid):
        if segid == segid2_name:
            resnum_list.append(float(readpdb_obj.atm_resseqnum[ind]))
    resnum_range = np.unique(np.asarray(resnum_list, dtype='float'))
    for resnum in resnum_range:
        xcoor = []
        ycoor = []
        zcoor = []
        for ind, segid in enumerate(readpdb_obj.atm_segid):
            if segid == segid2_name:
                if resnum == float(readpdb_obj.atm_resseqnum[ind]):
                    xcoor.append(readpdb_obj.atm_xcoor[ind])
                    ycoor.append(readpdb_obj.atm_ycoor[ind])
                    zcoor.append(readpdb_obj.atm_zcoor[ind])
        segid2_com = centerofmass(xcoor, ycoor, zcoor)
        distance = distance_between_two_points(segid1_com, segid2_com)
        distance_list.append(distance)


    # for resnum in resnum_range:
    #     xcoor = []
    #     ycoor = []
    #     zcoor = []
    #     for ind, segid in enumerate(readpdb_obj.atm_segid):
    #         if segid == segid2_name:
    #             if resnum == float(readpdb_obj.atm_resseqnum[ind]):
    #                 xcoor.append(readpdb_obj.atm_xcoor[ind])
    #                 ycoor.append(readpdb_obj.atm_ycoor[ind])
    #                 zcoor.append(readpdb_obj.atm_zcoor[ind])
    #     segid2_com = centerofmass(xcoor, ycoor, zcoor)
    #     distance = distance_between_two_points(segid1_com, segid2_com)
    #     distance_list.append(distance)


    return distance_list, resnum_range

def plot_distance_histogram_1(distance_list, dirpath):
    fig, ax1 = plt.subplots()

    ax1.hist(distance_list, bins='auto', color='lightblue')
    ax1.set_xlabel('Centre of mass distance (Angstrom)')
    ax1.set_ylabel('Count')

    ax2 = ax1.twinx()
    ax2.hist(distance_list, bins='auto', color='red', histtype='step', cumulative=True)
    ax2.set_ylabel('C.F')

    plt.savefig(os.path.join(dirpath, 'COM_Distance_plot.pdf'))
    plt.close()

def generate_pdb_models_with_distance_cutoff(dist_cut_off, segid1_name_list, segid_name, dist_list, res_num, readpdb_object, dirpath):
    list_of_fields_distcutoff = []
    del_res = []
    for seg in segid1_name_list:
        for ind_1, segid_1 in enumerate(readpdb_object.atm_segid):
            if segid_1 == seg:
                list_of_fields_distcutoff.append(
                    [readpdb_object.atm[ind_1], readpdb_object.atm_sernum[ind_1], readpdb_object.atm_name[ind_1],
                     readpdb_object.atm_altloc[ind_1],
                     readpdb_object.atm_resname[ind_1],
                     readpdb_object.atm_chainid[ind_1], readpdb_object.atm_resseqnum[ind_1],
                     readpdb_object.atm_cod[ind_1],
                     readpdb_object.atm_xcoor[ind_1], readpdb_object.atm_ycoor[ind_1],
                     readpdb_object.atm_zcoor[ind_1], readpdb_object.atm_occ[ind_1],
                     readpdb_object.atm_tempf[ind_1], segid_1, readpdb_object.atm_elem[ind_1],
                     readpdb_object.atm_charge[ind_1]])
    for index, (res, dist) in enumerate(zip(res_num, dist_list)):
        for ind, segid in enumerate(readpdb_object.atm_segid):
            if segid == segid_name:
                if dist <= dist_cut_off:
                    if float(readpdb_object.atm_resseqnum[ind]) == res:
                        list_of_fields_distcutoff.append(
                            [readpdb_object.atm[ind], readpdb_object.atm_sernum[ind], readpdb_object.atm_name[ind],
                             readpdb_object.atm_altloc[ind],
                             readpdb_object.atm_resname[ind],
                             readpdb_object.atm_chainid[ind], readpdb_object.atm_resseqnum[ind],
                             readpdb_object.atm_cod[ind],
                             readpdb_object.atm_xcoor[ind], readpdb_object.atm_ycoor[ind],
                             readpdb_object.atm_zcoor[ind], readpdb_object.atm_occ[ind],
                             readpdb_object.atm_tempf[ind], segid, readpdb_object.atm_elem[ind],
                             readpdb_object.atm_charge[ind]])
                if dist >= dist_cut_off:
                    if float(readpdb_object.atm_resseqnum[ind]) == res:
                        del_res.append(res)

    generate_del_res_stream_file('MEMB', del_res, dirpath)

    write_pdb(list_of_fields_distcutoff, str(dist_cut_off)+'_angstrom_cutoff_model.pdb', dirpath, elem_charge=True)

def generate_del_segid_res_stream_file(segid_list, resnum_list, dirpath, fname=None):
    if fname == None:
        fname = 'del_res.str'
    uniq_resnum_list, unique_index = np.unique(resnum_list, return_index=True)
    uniq_segid_list = np.array(segid_list)[unique_index]
    with open(os.path.join(dirpath, fname), 'w') as delres:
        for ind, (uniq_segid, uniq_resnum) in enumerate(zip(uniq_segid_list, uniq_resnum_list)):
            line = 'DELETe ATOM SELE SEGId {} .and. RESId {} end\n'.format(uniq_segid, str(int(uniq_resnum)))
            delres.write(line)


def gen_del_segid_stream_file(segid_list, dirpath, fname=None):
    if fname == None:
        fname  = 'del_res.str'
    uniq_segid_list = np.unique(segid_list)
    with open(os.path.join(dirpath, fname), 'w') as delres:
        for ind, uniq_segid in enumerate(uniq_segid_list):
            line = 'DELETe ATOM SELE SEGId {} end\n'.format(uniq_segid)
            delres.write(line)
        delres.close()


def generate_del_res_stream_file(segid, resnum_list, dirpath, fname=None):
    if fname == None:
        fname = 'del_res.str'
    uniq_resnum_list = np.unique(resnum_list)
    with open(os.path.join(dirpath, fname), 'w') as delres:
        for ind in range(len(uniq_resnum_list)):
            line = 'DELETe ATOM SELE SEGId {} .and. RESId {} end\n'.format(segid, str(int(uniq_resnum_list[ind])))
            delres.write(line)


def generate_crd_models_with_distance_cutoff(dist_cut_off, segid1_name_list, segid_name, dist_list, res_num, readcrd_object, dirpath):
    list_of_fields_distcutoff = []
    del_res = []
    for seg in segid1_name_list:
        for ind_1, segid_1 in enumerate(readcrd_object.atm_segid):
            if segid_1 == seg:
                list_of_fields_distcutoff.append(
                    [readcrd_object.atm_sernum[ind_1], readcrd_object.resnum_relative[ind_1],
                     readcrd_object.atm_resname[ind_1], readcrd_object.atm_name[ind_1],
                     readcrd_object.atm_xcoor[ind_1], readcrd_object.atm_ycoor[ind_1],
                     readcrd_object.atm_zcoor[ind_1], segid_1, readcrd_object.atm_resseqnum[ind_1],
                     readcrd_object.atm_charge[ind_1]])
    for index, (res, dist) in enumerate(zip(res_num, dist_list)):
        for ind, segid in enumerate(readcrd_object.atm_segid):
            if segid == segid_name:
                if dist <= dist_cut_off:
                    if float(readcrd_object.atm_resseqnum[ind]) == res:
                        list_of_fields_distcutoff.append(
                            [readcrd_object.atm_sernum[ind], readcrd_object.resnum_relative[ind],
                             readcrd_object.atm_resname[ind], readcrd_object.atm_name[ind],
                             readcrd_object.atm_xcoor[ind], readcrd_object.atm_ycoor[ind],
                             readcrd_object.atm_zcoor[ind], segid, readcrd_object.atm_resseqnum[ind],
                             readcrd_object.atm_charge[ind]])
                if dist >= dist_cut_off:
                    if float(readcrd_object.atm_resseqnum[ind]) == res:
                        del_res.append(res)

    generate_del_res_stream_file('MEMB', del_res, dirpath)

    write_crd(list_of_fields_distcutoff, str(dist_cut_off)+'_angstrom_cutoff_model.crd', dirpath)


# def generate_writepdb_list_of_fields_with_combinations_of_resseqnum(res_num, number):
#     combination = combinations(res_num, number)

def make_new_dir(dirpath, foldername):
    new_path = os.path.join(dirpath, foldername)
    if not os.path.isdir(new_path): os.makedirs(new_path)
    return new_path


def get_files(dirpath, endid='.pdb'):
    """
    get files from a directory ending with endid string
    :param dirpath: directory path
    :param endid: str
    :return: list of files
    """
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith(endid)]
    return files


def generate_models_with_number_of_lipids_plus_distance(number, distance, segid1_name_list, segid_name, dist_list, res_num, readpdb_object, dirpath):
    res_above_dist_cutoff, res_below_dist_cutoff = [], []
    dist_list_above_dist_cutoff, dist_list_below_dist_cutoff = [], []
    for index, (res, dist) in enumerate(zip(res_num, dist_list)):
        if dist >= distance:
            res_above_dist_cutoff.append(res)
            dist_list_above_dist_cutoff.append(dist)
        else:
            res_below_dist_cutoff.append(res)
            dist_list_below_dist_cutoff.append(dist)

    comb_lipids = combinations(res_above_dist_cutoff, number)

    num_lipids_base = len(res_below_dist_cutoff)
    total_num_lipids = len(res_below_dist_cutoff) + number

    out_folder_name = str(total_num_lipids)+'_numlipids'
    output_path = make_new_dir(dirpath, out_folder_name)


    comb_count = 0
    for combo in comb_lipids:
        print(combo)
        list_of_fields = []
        del_res = []
        for seg in segid1_name_list:
            for ind_1, segid_1 in enumerate(readpdb_object.atm_segid):
                if segid_1 == seg:
                    list_of_fields.append(
                        [readpdb_object.atm[ind_1], readpdb_object.atm_sernum[ind_1], readpdb_object.atm_name[ind_1],
                         readpdb_object.atm_altloc[ind_1],
                         readpdb_object.atm_resname[ind_1],
                         readpdb_object.atm_chainid[ind_1], readpdb_object.atm_resseqnum[ind_1],
                         readpdb_object.atm_cod[ind_1],
                         readpdb_object.atm_xcoor[ind_1], readpdb_object.atm_ycoor[ind_1],
                         readpdb_object.atm_zcoor[ind_1], readpdb_object.atm_occ[ind_1],
                         readpdb_object.atm_tempf[ind_1], segid_1, readpdb_object.atm_elem[ind_1],
                         readpdb_object.atm_charge[ind_1]])

        if res_below_dist_cutoff:
           for num, res_ in enumerate(res_below_dist_cutoff):
               for ind_base, segid_base in enumerate(readpdb_object.atm_segid):
                   if segid_base == segid_name:
                       if float(readpdb_object.atm_resseqnum[ind_base]) == res_:
                           list_of_fields.append([readpdb_object.atm[ind_base], readpdb_object.atm_sernum[ind_base], readpdb_object.atm_name[ind_base],
                         readpdb_object.atm_altloc[ind_base],
                         readpdb_object.atm_resname[ind_base],
                         readpdb_object.atm_chainid[ind_base], readpdb_object.atm_resseqnum[ind_base],
                         readpdb_object.atm_cod[ind_base],
                         readpdb_object.atm_xcoor[ind_base], readpdb_object.atm_ycoor[ind_base],
                         readpdb_object.atm_zcoor[ind_base], readpdb_object.atm_occ[ind_base],
                         readpdb_object.atm_tempf[ind_base], segid_base, readpdb_object.atm_elem[ind_base],
                         readpdb_object.atm_charge[ind_base]])



        for index, resnum in enumerate(combo):
            for ind, segid in enumerate(readpdb_object.atm_segid):
                if segid == segid_name:
                    if float(readpdb_object.atm_resseqnum[ind]) == resnum:
                        list_of_fields.append([readpdb_object.atm[ind], readpdb_object.atm_sernum[ind], readpdb_object.atm_name[ind],
                             readpdb_object.atm_altloc[ind],
                             readpdb_object.atm_resname[ind],
                             readpdb_object.atm_chainid[ind], readpdb_object.atm_resseqnum[ind],
                             readpdb_object.atm_cod[ind],
                             readpdb_object.atm_xcoor[ind], readpdb_object.atm_ycoor[ind],
                             readpdb_object.atm_zcoor[ind], readpdb_object.atm_occ[ind],
                             readpdb_object.atm_tempf[ind], segid, readpdb_object.atm_elem[ind],
                             readpdb_object.atm_charge[ind]])
                    if float(readpdb_object.atm_resseqnum[ind]) != resnum:
                        del_res.append(float(readpdb_object.atm_resseqnum[ind]))

        write_pdb(list_of_fields, '_numlipids_'+str(total_num_lipids)+'_combination_'+str(comb_count)+'_distance_'+
                  str(distance)+'.pdb', output_path, elem_charge=True)
        generate_del_res_stream_file('MEMB', del_res, output_path, 'numlipids_'+str(total_num_lipids)+'_comb_'+str(comb_count)+'_distance_'+str(distance)+'.str')
        comb_count += 1

def generate_models_with_increasing_distance_num_lipids(num_lipids, segid1_name_list, segid_name, dist_list, res_num, readpdb_object, dirpath):

    sort_index = np.argsort(dist_list)
    dist_list_sort = np.array(dist_list)[sort_index]
    res_num_sort = np.array(res_num)[sort_index]
    list_of_fields_distcutoff = []
    del_res = []
    for seg in segid1_name_list:
        for ind_1, segid_1 in enumerate(readpdb_object.atm_segid):
            if segid_1 == seg:
                list_of_fields_distcutoff.append(
                    [readpdb_object.atm[ind_1], readpdb_object.atm_sernum[ind_1], readpdb_object.atm_name[ind_1],
                     readpdb_object.atm_altloc[ind_1],
                     readpdb_object.atm_resname[ind_1],
                     readpdb_object.atm_chainid[ind_1], readpdb_object.atm_resseqnum[ind_1],
                     readpdb_object.atm_cod[ind_1],
                     readpdb_object.atm_xcoor[ind_1], readpdb_object.atm_ycoor[ind_1],
                     readpdb_object.atm_zcoor[ind_1], readpdb_object.atm_occ[ind_1],
                     readpdb_object.atm_tempf[ind_1], segid_1, readpdb_object.atm_elem[ind_1],
                     readpdb_object.atm_charge[ind_1]])
    for index in range(num_lipids):
        for ind, segid in enumerate(readpdb_object.atm_segid):
            if segid == segid_name:
                if float(readpdb_object.atm_resseqnum[ind]) == res_num_sort[index]:
                    list_of_fields_distcutoff.append(
                            [readpdb_object.atm[ind], readpdb_object.atm_sernum[ind], readpdb_object.atm_name[ind],
                             readpdb_object.atm_altloc[ind],
                             readpdb_object.atm_resname[ind],
                             readpdb_object.atm_chainid[ind], readpdb_object.atm_resseqnum[ind],
                             readpdb_object.atm_cod[ind],
                             readpdb_object.atm_xcoor[ind], readpdb_object.atm_ycoor[ind],
                             readpdb_object.atm_zcoor[ind], readpdb_object.atm_occ[ind],
                             readpdb_object.atm_tempf[ind], segid, readpdb_object.atm_elem[ind],
                             readpdb_object.atm_charge[ind]])

    for index in range(num_lipids, len(res_num_sort)):
        del_res.append(res_num_sort[index])

    generate_del_res_stream_file('MEMB', del_res, dirpath, fname='del_res_'+str(num_lipids)+'.str')

    write_pdb(list_of_fields_distcutoff, str(num_lipids) + '_incr_dist_model.pdb', dirpath, elem_charge=True)



def generate_models_with_delete_resseqnum(readpdb_object, resseqnum, dirpath):
    list_of_fields = []
    for ind, res_num in enumerate(readpdb_object.atm_resseqnum):
        if float(res_num) != resseqnum:
            list_of_fields.append([readpdb_object.atm[ind], readpdb_object.atm_sernum[ind], readpdb_object.atm_name[ind],
                             readpdb_object.atm_altloc[ind],
                             readpdb_object.atm_resname[ind],
                             readpdb_object.atm_chainid[ind], res_num,
                             readpdb_object.atm_cod[ind],
                             readpdb_object.atm_xcoor[ind], readpdb_object.atm_ycoor[ind],
                             readpdb_object.atm_zcoor[ind], readpdb_object.atm_occ[ind],
                             readpdb_object.atm_tempf[ind], readpdb_object.atm_segid[ind], readpdb_object.atm_elem[ind],
                             readpdb_object.atm_charge[ind]])
    original_pdb_file_name = readpdb_object.pdbfile.split('.')[0]
    pdb_fname = original_pdb_file_name+'_del_.pdb'
    write_pdb(list_of_fields, pdb_fname, dirpath, elem_charge=True)



def generate_pdb_with_select_segids(dirpath_with_pdbs, list_of_segids, output_dir):
    """
    generate pdb with including only select segids and also write stream file
    :param list_of_pdb_fpaths: list of pdb fpaths
    :param list_of_segids: list of segids
    :param output_dir: where to output files
    :return: none. writes files -> pdb and stream files
    """

    pdb_files = get_files(dirpath_with_pdbs, endid='.pdb')

    output_path = make_new_dir(dirpath_with_pdbs, output_dir)

    for index, pdb_file in enumerate(pdb_files):

        print ('Generating '+str(index+1)+' of '+str(len(pdb_files))+' files')

        readpdb_object = ParsePDB(pdb_file, dirpath_with_pdbs).read()

        list_of_fields = []


        for seg in list_of_segids:
            for ind_1, segid_1 in enumerate(readpdb_object.atm_segid):
                if segid_1 == seg:
                    list_of_fields.append(
                        [readpdb_object.atm[ind_1], readpdb_object.atm_sernum[ind_1], readpdb_object.atm_name[ind_1],
                         readpdb_object.atm_altloc[ind_1],
                         readpdb_object.atm_resname[ind_1],
                         readpdb_object.atm_chainid[ind_1], readpdb_object.atm_resseqnum[ind_1],
                         readpdb_object.atm_cod[ind_1],
                         readpdb_object.atm_xcoor[ind_1], readpdb_object.atm_ycoor[ind_1],
                         readpdb_object.atm_zcoor[ind_1], readpdb_object.atm_occ[ind_1],
                         readpdb_object.atm_tempf[ind_1], segid_1, readpdb_object.atm_elem[ind_1],
                         readpdb_object.atm_charge[ind_1]])

        # write new pdb file
        original_pdb_file_name = readpdb_object.pdbfile.split('.')[0]
        pdb_fname = original_pdb_file_name + '_del_.pdb'
        write_pdb(list_of_fields, pdb_fname, output_path, elem_charge=True)

        # write del res stream file
        # gen_del_segid_stream_file(list_of_del_segids, output_path, original_pdb_file_name+'_del_res.str')


if __name__=='__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\charmm_gui_iapp_nanodisc\iappdimer_peripheral_2\19"
    pdbfile = "step7_4.pdb"

    mod = ParsePDB(pdbfile, dirpath).read()
    iappdimer_com = compute_com_segid(mod, ['PROA', 'PROB'])
    dist_list, res_num = compute_distance_list(mod, iappdimer_com, 'MEMB')

    # plt.hist(dist_list, bins='auto')
    # plt.savefig(os.path.join(dirpath, 'dist_list.png'))
    # plt.close()

    # generate_models_with_number_of_lipids_plus_distance(1, 10, ['PROA', 'PROB'], 'MEMB', dist_list, res_num, mod, dirpath)

    #
    # generate_pdb_models_with_distance_cutoff(15, ['PROA', 'PROB'], 'MEMB', dist_list, res_num, mod, dirpath)

    generate_models_with_increasing_distance_num_lipids(10, ['PROA', 'PROB'], 'MEMB', dist_list, res_num, mod, dirpath)

    # numlipids11_pdbfile = '11_incr_dist_model.pdb'
    # numlipids11_pdb = ParsePDB(numlipids11_pdbfile, dirpath).read()
    # iappdimer_com = compute_com_segid(numlipids11_pdb, ['PROA', 'PROB'])
    # dist_list, res_num = compute_distance_list(numlipids11_pdb, iappdimer_com, 'MEMB')
    #
    # generate_models_with_number_of_lipids_plus_distance(5, 5, ['PROA', 'PROB'], 'MEMB', dist_list, res_num, numlipids11_pdb, dirpath)

    # dirpath_with_pdbs = r"C:\Users\sugyan\Documents\MembraneMDfiles\charmm_gui_iapp_nanodisc\iappdimer_integral\19\step7_structs"
    # generate_pdb_with_select_segids(dirpath_with_pdbs, ['PROA', 'PROB'], 'output')
