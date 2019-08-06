# calculate dipole moment for each cluster based on the start and end residue list

import os
import pickle
import pdb_utility
from pdb_utility import ParsePDB, calculate_helix_dipole#, centerofmass, distance_between_two_points

def get_files(dirpath, endid, startid=None):

    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith(endid)]
    if startid:
        files = [x for x in files if x.startswith(startid)]
    return files


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


class DipoleMoment(object):
    """
    object to store dipole moment information from pdb files
    """

    def __init__(self):

        self.res_start_end_list = []
        self.pdb_fpath = []
        self.dipole_moment_dict = []


def write_dipole_moment_output(dipole_moment_object, dirpath):
    """
    with the dipole moment object file write the essential output
    :param dipole_moment_object: dipole moment object
    :return: void. writes the output file
    """

    header = 'pdb_fpath,start_res_num,end_res_num,dipole,dipole_per_res\n'

    data_string = ''

    for ind, (pdb_file, dipole_dict) in enumerate(zip(dipole_moment_object.pdb_fpath, dipole_moment_object.dipole_moment_dict)):

        line = '{},{},{},{},{}\n'.format(pdb_file,
                                         dipole_dict['start_res_num'],
                                         dipole_dict['end_res_num'],
                                         dipole_dict['dipole'],
                                         dipole_dict['dipole_per_res'])
        data_string += line


    output_string = header + data_string

    outfname = os.path.join(dirpath, 'dipole_moment_list.csv')

    with open(outfname, 'w') as outfile:
        outfile.write(output_string)
        outfile.close()


def gen_dipole_moment_obj(dirpath, pdb_file_list, res_start_end_list):
    """
    read dssp list obj
    :param res_start_end: list [start_res, end_res]
    :param pdb_fpath_list: list of pdb fpaths for which to calculate the dipole moment
    :param dirpath: directory to save object file
    :return: dipole_moment_object
    """

    dipole_moment_obj = DipoleMoment()

    for ind, res_start_end in enumerate(res_start_end_list):
        dipole_moment_obj.res_start_end_list.append(res_start_end)
        for index, pdb_file in enumerate(pdb_file_list):
            pdb_read = ParsePDB(pdb_file, dirpath).read()
            dipole_mom_dict = calculate_helix_dipole(res_start_end[0], res_start_end[1], pdb_read)
            dipole_moment_obj.pdb_fpath.append(pdb_file)
            dipole_moment_obj.dipole_moment_dict.append(dipole_mom_dict)


    write_dipole_moment_output(dipole_moment_obj, dirpath)

    outname = os.path.join(dirpath, 'dipole_moment_object.obj')

    save_object_to_pickle(dipole_moment_obj, outname)

    return dipole_moment_obj


if __name__ == '__main__':

    dirpath = os.getcwd()
    pdb_file_list = get_files(dirpath, endid='.pdb', startid='repid')

    start_end_res_list = [[19, 20], [13, 18], [30, 36], [21, 27], [38, 43], [44, 49], [53, 58]]

    dipole_moment_calc = gen_dipole_moment_obj(dirpath, pdb_file_list, start_end_res_list)
