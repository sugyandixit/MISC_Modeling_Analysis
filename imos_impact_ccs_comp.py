# match ccs from pdb file for imos and impact
# do a correlation (linear or other functions)
# save outputs and plot

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def plot_ccs_comp(imp_ccs_arr, imos_ccs_arr, linreg):
    """
    plot imp vs imos ccs with linreg line and text
    :param imp_ccs_arr: imp_ccs_arr
    :param imos_ccs_arr: imos_ccs_arr
    :param linreg: linreg
    :return: void. outfile
    """

    x_lin_fit = np.linspace(np.min(imp_ccs_arr), np.max(imp_ccs_arr), 500)
    y_lin_fit = x_lin_fit * linreg[0] + linreg[1]

    plt.scatter(imp_ccs_arr, imos_ccs_arr, marker='o', color='black')
    plt.plot(x_lin_fit, y_lin_fit, ls='--', color='red')
    plt.xlabel('IMPACT_CCS_TM')
    plt.ylabel('IMOS_CCS_TM')


def write_ccs_comp(clust_num, pdb_fpath, imp_ccs_arr, imos_ccs_arr, linreg, rmse, output_dir, outid):

    linreg_line0 = '#xdata,imp_ccs_arr\n#ydata,imos_ccs_arr\n'
    linreg_line1 = '#slope,'+str(linreg[0]) + '\n'
    linreg_line2 = '#intercept,' + str(linreg[1]) + '\n'
    linreg_line3 = '#r2,'+str(linreg[2]**2) + '\n'
    linreg_line4 = '#std_err,'+str(linreg[4]) + '\n'

    rmse_line = '#RMSE(%),'+str(rmse)+'\n'

    linreg_line = linreg_line0 + linreg_line1 + linreg_line2 + linreg_line3 + linreg_line4

    header = 'clust_num,pdb_fpath,imp_ccs_tm,imos_ccs_tm\n'

    data_string = ''

    for num in range(len(imos_ccs_arr)):
        line = '{},{},{},{}\n'.format(clust_num[num],
                                      pdb_fpath[num],
                                      imp_ccs_arr[num],
                                      imos_ccs_arr[num])
        data_string += line

    output_string = linreg_line + rmse_line + header + data_string

    with open(os.path.join(output_dir, 'imp_imos_ccs_'+outid+'.csv'), 'w') as outfile:
        outfile.write(output_string)
        outfile.close()


def rms_error(predict_values, true_values):
    """
    Calculate the percent error and computes the rms error in %
    :param predict_values: predict values
    :param true_values: true values
    :return: % rmse
    """
    diff = np.subtract(predict_values, true_values)
    percent_error = np.multiply(np.divide(diff, true_values), 100)
    percent_error_sq = np.square(percent_error)
    rmse = np.sqrt(np.divide(np.sum(percent_error_sq), len(percent_error_sq)))
    return rmse




def sort_ccs_imos_impact(list_of_impact_imos_file_pair, output_dir, outid, image_type='.pdf'):
    """
    sort ccs based on same pdb file name
    :param list_of_impact_imos_file_pair: list of impact imos file pair
    :return: void. Ouptut file for matched ccs, correlation, and plot imos vs impact ccs
    """

    clust_num = []
    imp_pdb_fpath = []
    imos_ccs_arr = []
    impact_ccs_arr = []

    for ind, file_pair_list in enumerate(list_of_impact_imos_file_pair):
        imos_file, impact_file = file_pair_list
        imos_df = pd.read_csv(imos_file)
        impact_df = pd.read_csv(impact_file)

        for index1, (imos_pdb_fpath, imos_ccs) in enumerate(zip(imos_df['file'].values, imos_df['CCS_TM'].values)):
            imos_pdb_fname = str(os.path.split(imos_pdb_fpath)[1])
            imos_pdb_fname = imos_pdb_fname.split('_out.pdb')[0] + '.pdb'
            for index2, (impact_pdb_fpath, impact_ccs) in enumerate(zip(impact_df['pdb_fname'].values, impact_df['CCS_TM'].values)):
                impact_pdb_fname = str(os.path.split(impact_pdb_fpath)[1])
                if imos_pdb_fname == impact_pdb_fname:
                    clust_num.append(impact_df['clust_num'].values[index2])
                    imp_pdb_fpath.append(impact_pdb_fpath)
                    imos_ccs_arr.append(imos_ccs)
                    impact_ccs_arr.append(impact_ccs)

    imos_ccs_arr = np.array(imos_ccs_arr)
    impact_ccs_arr = np.array(impact_ccs_arr)

    linreg = linregress(impact_ccs_arr, imos_ccs_arr)
    y_pred = impact_ccs_arr * linreg[0] + linreg[1]
    rmse = rms_error(y_pred, impact_ccs_arr)

    write_ccs_comp(clust_num, imp_pdb_fpath, impact_ccs_arr, imos_ccs_arr, linreg, rmse, output_dir, outid)

    plot_ccs_comp(impact_ccs_arr, imos_ccs_arr, linreg)
    plt.savefig(os.path.join(output_dir, 'imp_imos_ccs_plot_'+ outid + image_type))
    plt.close()

    print('heho')


if __name__ == '__main__':

    p1_imos_file = r"C:\Users\sugyan\Documents\IMOS\savefolder\serf_mod_structs_random_selection\IMOS_CCS_out.csv"
    p1_imp_file = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\rmsd_mat\cluster_manual_pdb_ccs_randomselection.csv"

    p2_imos_file = r"C:\Users\sugyan\Documents\IMOS\savefolder\serf_nomod_structs_random_selection\IMOS_CCS_out.csv"
    p2_imp_file = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_nomod\rmsdmat\cluster_manual_pdb_ccs_randomselection.csv"

    output_dir = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\_Exp_CCS"

    imp_imos_list = [[p2_imos_file, p2_imp_file]]#,[p2_imos_file, p2_imp_file]]

    sort_ccs_imos_impact(imp_imos_list, output_dir, outid='nomod')