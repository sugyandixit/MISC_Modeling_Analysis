import os
import pickle
import itertools
import numpy as np
from sklearn.metrics import jaccard_similarity_score
from scipy.spatial.distance import jaccard, squareform, pdist
import scipy.cluster.hierarchy as hclust
import matplotlib.pyplot as plt
from generate_contact_map_list_object import ContactMap
import time


def write_jac_dist(pdb_comb_list, jac_dist_array, dist_cut_off, dirpath):
    outlines = ''
    header='#struct1,struct2,jac_dist\n'
    outlines += header

    for index, (pdb_comb_, jac_dist) in enumerate(zip(pdb_comb_list, jac_dist_array)):
        line = '{},{},{}\n'.format(pdb_comb_[0], pdb_comb_[1], jac_dist)
        outlines += line

    with open(os.path.join(dirpath, 'jac_dist_ct_map_'+ str(dist_cut_off)+'.csv'), 'w') as outfile:
        outfile.write(outlines)
        outfile.close()


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
    ct_map_fname = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\contact_mat_10.obj"
    dirpath = os.path.split(ct_map_fname)[0]
    ct_map_obj = load_pickle_object(ct_map_fname)
    jac_score_ = []
    jac_dist_ = []
    jac_sim_ = []
    jac_dist_sim_ = []
    pdb_comb_list = []
    jac_dist_scipy_ = []
    start_time1 = time.perf_counter()
    start_time2 = time.process_time()

    with open(os.path.join(dirpath, 'jac_dist_ct_map_' + str(ct_map_obj.dist_cut_off) + '.csv'), 'w') as outfile:
        header = '#struct1,struct2,jac_dist\n'
        outfile.write(header)

        for ind, (ct_comb, pdb_comb) in enumerate(zip(itertools.combinations(ct_map_obj.contact_mat, 2),
                                                      itertools.combinations(ct_map_obj.pdbname, 2))):

            print('for contact dist cut off ' + str(ct_map_obj.dist_cut_off) + ' angstrom calculating jaccard distance between ' + pdb_comb[0] + ' and ' + pdb_comb[1] + '    ...')

            ct_comb0_flat = ct_comb[0].flatten()
            ct_comb1_flat = ct_comb[1].flatten()

            try:
                jac_dist_scipy = jaccard(ct_comb0_flat, ct_comb1_flat)
                comb_files = pdb_comb
            except RuntimeWarning:
                print('invalid value encountered ... discarding the combination ...')
            # jac_dist_scipy_.append(jac_dist_scipy)
            # pdb_comb_list.append(comb_files)

            line = '{},{},{}\n'.format(comb_files[0], comb_files[1], jac_dist_scipy)
            outfile.write(line)


            elapsed_time = time.perf_counter() - start_time1
            cpu_process_time = time.process_time() - start_time2
            print('elapsed time ... ' + str(elapsed_time))
            print('cpu time ... ' + str(cpu_process_time))
            start_time1 = time.perf_counter()
            start_time2 = time.process_time()

        outfile.close()

        # jac_dist_scipy_ = np.array(jac_dist_scipy_)
        #
        # print('writing output ...')
        # write_jac_dist(pdb_comb_list, jac_dist_scipy_, ct_map_obj.dist_cut_off, dirpath)
        #
        #
        # print('heho')

        # os.path.join(os.getcwd(), 'pyrosetta_output/contact_mat_10.obj')