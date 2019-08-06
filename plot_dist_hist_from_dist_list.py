
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def get_files(dirpath, endid, startid=None):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith(endid)]
    if startid:
        files = [x for x in files if x.startswith(startid)]

    return files


def write_dist_output(file, dirpath, outid):

    header = 'res_num,dist_label,dist_avg,dist_std\n'

    data_string = ''

    df = pd.read_csv(os.path.join(dirpath, file))

    uniq_ref_res = np.unique(df['ref_res_num'].values)

    for ind, uniq_res in enumerate(uniq_ref_res):

        df_uniq_res = df[(df['ref_res_num'] == uniq_res)]

        df_dist = df_uniq_res.iloc[:, 5:]
        df_dist_values = df_dist.values
        df_dist_columns = df_dist.columns.values

        df_dist_t = df_dist_values.T

        for index, (arr, label) in enumerate(zip(df_dist_t, df_dist_columns)):
            mean_dist = np.mean(arr)
            std_dist = np.std(arr)

            line = '{},{},{},{}\n'.format(uniq_res,
                                          label,
                                          mean_dist,
                                          std_dist)
            data_string += line

    output_string = header + data_string

    outname = outid + 'dist_stat_out.csv'

    with open(os.path.join(dirpath, outname), 'w') as outfile:
        outfile.write(output_string)
        outfile.close()



def plot_dist_hist(file, dirpath, outid):

    df = pd.read_csv(os.path.join(dirpath, file))

    uniq_ref_res = np.unique(df['ref_res_num'].values)

    outname = outid + 'dist_hist.pdf'

    with PdfPages(os.path.join(dirpath, outname)) as pdf:

        for ind, uniq_res in enumerate(uniq_ref_res):

            df_uniq_res = df[(df['ref_res_num'] == uniq_res)]

            df_hist = df_uniq_res.iloc[:, 5:]
            df_hist_values = df_hist.values
            df_hist_columns = df_hist.columns.values

            df_hist_t = df_hist_values.T

            for index, (arr, label) in enumerate(zip(df_hist_t, df_hist_columns)):

                plt.hist(arr, bins='auto', density=False, label=label, alpha=1)
                plt.title('RES_NUM_'+str(uniq_res))

            plt.legend(loc='best', fontsize='small')
            plt.xlim(left=0.5, right=25)
            plt.xlabel('dist')
            plt.ylabel('Count')
            pdf.savefig()
            plt.close()


if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\rmsd_mat"

    files = get_files(dirpath, endid='dist_list_output.csv')

    for ind, file in enumerate(files):

        outid = str(file).split('.csv')[0]

        write_dist_output(file, dirpath, outid)
        plot_dist_hist(file, dirpath, outid)