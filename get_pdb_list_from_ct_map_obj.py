# list_of_pdb_files_from_contact_map
import os
import pickle
import numpy as np
import pandas as pd
from generate_contact_map_list_object import ContactMap


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


def get_contact_obj_file(dirpath, ct_dist=5):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('.obj')]
    files = [x for x in files if x.startswith('contact_mat_' + str(ct_dist))]
    return files


def write_output_file(ct_dist, pdb_list, dir_list, pyrosetta_score, dirpath):
    dataframe = pd.DataFrame({'dirpath':dir_list, 'pdb_name':pdb_list, 'rosetta_score':pyrosetta_score})
    dataframe_unique = dataframe.drop_duplicates(subset='pdb_name')
    dataframe_unique.to_csv(os.path.join(dirpath, 'pdb_list_ct_map_'+str(ct_dist)+'.csv'), sep=',', index=False)
    print('heho')


if __name__ == '__main__':
    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\cluster_1"
    # dirpath = os.path.join(os.getcwd(), 'pyrosetta_output')
    # pyros_score_out_file = os.path.join(os.getcwd(), 'pyrosetta_filter_output.csv')
    pyros_score_out_file = os.path.join(dirpath, 'pyrosetta_filter_output.csv')
    pyros_df = pd.read_csv(pyros_score_out_file, sep=',')
    ct_dist = 5
    contact_obj_files = get_contact_obj_file(dirpath, ct_dist)
    for ind, file in enumerate(contact_obj_files):
        obj_file = load_pickle_object(os.path.join(dirpath, file))
        pdb_list = obj_file.pdbname
        dirlist = obj_file.dirpath
        pdb_fname_pyros = pyros_df['#fname']
        ros_score = []
        for ind, pdbname in enumerate(pdb_list):
            score = pyros_df[pyros_df['#fname'] == pdbname]['total_score'].values[0]
            ros_score.append(score)

        write_output_file(ct_dist, pdb_list, dirlist, ros_score, dirpath)