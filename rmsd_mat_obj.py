import os
import pickle
import time
import pandas as pd
import numpy as np
import itertools
import rmsd_calculate
from rmsd_calculate import calculate_rmsd


class RMSDMAT(object):
    """
    Object to hold RMSD matrix, input pdb files, and combinations
    """
    def __init__(self):
        """
        initialize the object
        """
        self.input_pdb_files = []
        self.rmsd_list = []

def save_object_to_pickle(obj, file_path):
    with open(file_path, 'wb') as out_obj:
        pickle.dump(obj, out_obj)
        out_obj.close()


def get_ct_map_pdb_list_file(dirpath, cdist):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('.csv')]
    files = [x for x in files if x.startswith('pdb_list_ct_map_'+str(cdist))]
    return files

def get_pdb_files_in_dir(dirpath):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('.pdb')]
    return files


if __name__ == '__main__':
    # dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\cluster_0"
    dirpath = os.getcwd()
    start_time_ = time.perf_counter()
    cpu_start_time_ = time.process_time()
    cdist = 5
    # pdb_list_file = get_ct_map_pdb_list_file(dirpath, cdist)


    rmsd_mat_obj = RMSDMAT()

    all_pdb_files = get_pdb_files_in_dir(dirpath)
    # sample_frame = np.linspace(9995, 20000, 2001, dtype=int)
    # for framenum in sample_frame:
    #     for file in all_pdb_files:
    #         file_chars = str(file).split('.pdb')[0]
    #         frame_num = int(file_chars.split('_')[-1])
    #         if frame_num == framenum:
    #             rmsd_mat_obj.input_pdb_files.append(os.path.join(dirpath, file))

    for file in all_pdb_files:
        rmsd_mat_obj.input_pdb_files.append(os.path.join(dirpath, file))

    print(rmsd_mat_obj.input_pdb_files)
    print(len(rmsd_mat_obj.input_pdb_files))


    start_time = time.perf_counter()
    cpu_start_time = time.process_time()

    for ind, (pdb1, pdb2) in enumerate(itertools.combinations(rmsd_mat_obj.input_pdb_files, 2)):
        print ('Calculating rmsd between ', pdb1, ' & ', pdb2, ' ... ')
        rmsd = calculate_rmsd(pdb1, pdb2)
        check1_et = time.perf_counter()
        check1_cput = time.process_time()
        print('elapsed_time = ', check1_et - start_time)
        print('process_time = ', check1_cput - cpu_start_time)
        start_time = check1_et
        cpu_start_time = check1_cput
        rmsd_mat_obj.rmsd_list.append(rmsd)
    rmsd_mat_obj.rmsd_list = np.array(rmsd_mat_obj.rmsd_list)
    rmsd_mat_obj_fname = os.path.join(dirpath, 'rmsd_mat.obj')
    print('Saving object file ', rmsd_mat_obj_fname, '... ')
    save_object_to_pickle(rmsd_mat_obj, rmsd_mat_obj_fname)
    print('heho')
    print(' - - - ')
    print(' |   | ')
    print('  * *  ')
    print('   |   ')
    print('   |   ')
    print('   V   ')
    print('D O N E')
    print('total elapsed time = ', time.perf_counter() - start_time_)
    print ('total processed time = ', time.process_time() - cpu_start_time_)

