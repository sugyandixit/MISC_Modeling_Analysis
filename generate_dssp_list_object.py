import os
import pickle

import read_dssp_file
from read_dssp_file import DSSPObject, create_dssp_obj


def save_object_to_pickle(obj, file_path):
    with open(file_path, 'wb') as out_obj:
        pickle.dump(obj, out_obj)
        out_obj.close()


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


class DSSPListObj(object):
    """
    store all dssp object in a single object
    """

    def __init__(self):

        self.dssp_object = []



def get_dssp_file(dirpath):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('.dssp')]
    return files


def save_dssp_list_obj(dirpath):

    dssp_files = get_dssp_file(dirpath)

    dssp_list_obj = DSSPListObj()

    for file in dssp_files:
        dssp_fpath = os.path.join(dirpath, file)
        dssp_object = create_dssp_obj(dssp_fpath)
        dssp_list_obj.dssp_object.append(dssp_object)


    obj_name = 'dssp_list_object.obj'
    obj_fpath = os.path.join(dirpath, obj_name)
    save_object_to_pickle(dssp_list_obj, obj_fpath)


if __name__ == '__main__':

    # dirpath = os.getcwd()
    # save_dssp_list_obj(dirpath)
    dssp_list_obj_fpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\rmsd_mat\dssp_list_object.obj"
    dssp_list_obj = load_pickle_object(dssp_list_obj_fpath)
    print('heho')

