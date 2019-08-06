import os
import numpy as np
import pandas as pd


def get_files(dirpath, endid):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith(endid)]
    return files


def get_dipole_moment_summary(fpath):
    """
    look at unique res start or rest end num and get dipole,dipole/res avg and std
    :param fpath: filepath
    :return: void. Write the file
    """

    dirpath, fname = os.path.split(fpath)

    header = 'start_res_num,end_res_num,dipole_avg,dipole_std,dipole_res_avg,dipole_res_std\n'

    data_string = ''

    df = pd.read_csv(fpath)
    uniq_start_res = np.unique(df['start_res_num'].values)


    for ind, uniq_res in enumerate(uniq_start_res):

        df_uniq = df[(df['start_res_num']==uniq_res)]
        start_res = df_uniq['start_res_num'].values[0]
        end_res = df_uniq['end_res_num'].values[0]
        dipole_mean = np.mean(df_uniq['dipole'].values)
        dipole_std = np.std(df_uniq['dipole'].values)
        dipole_res_mean = np.mean(df_uniq['dipole_per_res'].values)
        dipole_res_std = np.std(df_uniq['dipole_per_res'].values)

        line = '{},{},{},{},{},{}\n'.format(start_res,
                                            end_res,
                                            dipole_mean,
                                            dipole_std,
                                            dipole_res_mean,
                                            dipole_res_std)
        data_string += line


    output_string = header + data_string

    outname = str(fname).split('.csv')[0] + '_summary.csv'

    with open(os.path.join(dirpath, outname), 'w') as outfile:
        outfile.write(output_string)
        outfile.close()

    print('heho')


if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\rmsd_mat"

    files = get_files(dirpath, endid='dipole_moment_list.csv')

    for file in files:

        fpath = os.path.join(dirpath, file)

        get_dipole_moment_summary(fpath)