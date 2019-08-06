import os
import pickle
import numpy as np
from rexmd_exchange_file_store_obj import RepExchOutObjCombined, RepExchOutObj


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

def get_exchange_ratio(rex_comb_obj, dirpath):
    """
    calculate the exchange ratio and write it out on a file
    :param rex_comb_obj:rex comb object file
    :param dirpath: directory to save the output file
    :return:void. Writes the file
    """

    mean_success_prob_temp = []
    std_success_prob_temp = []
    temp_ = []

    for ind, rex_obj in enumerate(rex_comb_obj.replica_exchange_obj):
        succes_prob = []
        for index, success_list in enumerate(rex_obj.success):
            if success_list == True:
                succes_prob.append(rex_obj.prob[index])

        mean_success_prob = np.mean(succes_prob)
        std_success_prob = np.std(succes_prob)

        mean_success_prob_temp.append(mean_success_prob)
        std_success_prob_temp.append(std_success_prob)

        temp_.append(rex_obj.temp_replica[0])

    header = 'temp,mean_success_prob,std_success_prob\n'

    data_string = ''

    for num in range(len(mean_success_prob_temp)):

        line = '{},{},{}\n'.format(temp_[num],
                                   mean_success_prob_temp[num],
                                   std_success_prob_temp[num])
        data_string += line


    output_string = header + data_string

    with open(os.path.join(dirpath, 'rex_exchange_prob.csv'), 'w') as outfile:
        outfile.write(output_string)
        outfile.close()

    print('heho')




if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\exch_files"
    rex_comb_obj_fpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\exch_files\rep_exch_comb_object.obj"
    rex_comb_obj = load_pickle_object(rex_comb_obj_fpath)
    get_exchange_ratio(rex_comb_obj, dirpath)
    print('heho')