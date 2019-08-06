import os
import pickle
import numpy as np


# store the replica exchange file in object

def save_object_to_pickle(obj, file_path):
    with open(file_path, 'wb') as out_obj:
        pickle.dump(obj, out_obj)
        out_obj.close()


def get_exch_files(dirpath):
    files = os.listdir(dirpath)
    files = [x for x in files if x.find('.exch_') >= 0]
    return files


class RepExchOutObjCombined(object):
    """
    stores a list of exchange out object
    """

    def __init__(self):
        self.replia_id = []
        self.replica_temps = None
        self.replica_exchange_obj = []


class RepExchOutObj(object):
    """
    create a replica exchange out object file
    """

    def __init__(self):
        """
        initialize the object with exchange file .exch
        """
        self.rep_id = None
        self.exchange_number = None
        self.step_number = None
        self.replica_num = None
        self.neighbor_num = None
        self.temp_replica = None
        self.temp_neighbor = None
        self.epot_replica = None
        self.epot_neighbor = None
        self.orig_tag = None
        self.new_tag = None
        self.prob = None
        self.rand = None
        self.tscale = None
        self.success = None


def assemble_rep_exch_comb_obj(dirpath):
    """
    ini
    :param dirpath:
    :param list_of_files:
    :return:
    """
    list_of_files = get_exch_files(dirpath)
    rep_exch_comb_obj = RepExchOutObjCombined()
    uniq_temps_each_replica = []

    for file in list_of_files:
        repx_exch_obj = read_exch_file(dirpath, file)
        rep_exch_comb_obj.replia_id.append(repx_exch_obj.rep_id)
        rep_exch_comb_obj.replica_exchange_obj.append(repx_exch_obj)
        repl_temp_uniq = np.unique(repx_exch_obj.temp_replica)
        uniq_temps_each_replica.append(repl_temp_uniq)

    rep_exch_comb_obj.replica_temps = np.unique(uniq_temps_each_replica)

    rep_exch_comb_obj_fname = 'rep_exch_comb_object.obj'
    save_object_to_pickle(rep_exch_comb_obj, os.path.join(dirpath, rep_exch_comb_obj_fname))


def read_exch_file(dirpath, fname):
    """
    read the file and store important information in lists
    :param fpath: file path
    :return: repexchoutobject
    """

    rep_id_num = int(fname.split('.exch_')[1])

    exch_num = []
    step_num = []
    replica_num = []
    temp_replica = []
    epot_replica = []
    neighbor_num = []
    temp_neighbor = []
    epot_neighbor = []
    orig_tag = []
    new_tag = []
    prob = []
    rand = []
    tscale = []
    success = []

    fpath = os.path.join(dirpath, fname)
    with open(fpath, 'r') as exchfile:
        exchfile_ = exchfile.readlines()
        for line in exchfile_:
            if line.startswith('REX>'):
                line = line.split('REX>')[1]
                line_chars = line.split()
                if line.startswith('EXCHANGE'):
                    exch_num.append(int(line_chars[2]))
                    step_num.append(int(line_chars[5]))
                if line.startswith('REPL'):
                    replica_num.append(int(line_chars[2]))
                    temp_replica.append(float(line_chars[5]))
                    epot_replica.append(float(line_chars[8]))
                if line.startswith('NEIGHBOR'):
                    neighbor_num.append(int(line_chars[2]))
                    temp_neighbor.append(float(line_chars[5]))
                    epot_neighbor.append(float(line_chars[8]))
                if line.startswith('ORIGINAL'):
                    orig_tag.append(int(line_chars[2]))
                    new_tag.append(int(line_chars[5]))
                if line.startswith('PROB'):
                    prob.append(float(line_chars[2]))
                    rand.append(float(line_chars[5]))
                    tscale.append(float(line_chars[8]))
                    if line_chars[11] == 'T':
                        success.append(True)
                    else:
                        success.append(False)

    rep_exch_obj = RepExchOutObj()
    rep_exch_obj.rep_id = rep_id_num
    rep_exch_obj.exchange_number = exch_num
    rep_exch_obj.step_number = step_num
    rep_exch_obj.replica_num = replica_num
    rep_exch_obj.temp_replica = temp_replica
    rep_exch_obj.epot_replica = epot_replica
    rep_exch_obj.neighbor_num = neighbor_num
    rep_exch_obj.temp_neighbor = temp_neighbor
    rep_exch_obj.epot_neighbor = epot_neighbor
    rep_exch_obj.orig_tag = orig_tag
    rep_exch_obj.new_tag = new_tag
    rep_exch_obj.prob = prob
    rep_exch_obj.rand = rand
    rep_exch_obj.tscale = tscale
    rep_exch_obj.success = success

    return rep_exch_obj


if __name__=='__main__':
    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis\exch_files"
    assemble_rep_exch_comb_obj(dirpath)