# cluster pdb match ccs

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt

def plot_ccs_mean_std(clust_num, mean_ccs, std_ccs, dirpath, outid):

    # x_data = [1 for i in range(len(clust_num))]

    # color_ = plt.get_cmap("tab10")
    color_ = ['red', 'royalblue', 'slategray', 'deepskyblue', 'darkslategray', 'green', 'gold', 'orange']

    for num in range(len(clust_num)):
        plt.errorbar(1, mean_ccs[num], yerr=std_ccs[num], marker='o', ls='none', capsize=2, capthick=1, ecolor='black', elinewidth=2,
                 markeredgecolor='black',
                 markersize=10,
                     markerfacecolor=color_[num],
                     label=str(clust_num[num]), alpha=1)
        plt.axhline(y=mean_ccs[num], ls='--', color=color_[num])

    plt.legend()
    plt.xlim((0,2))
    plt.ylim((800, 2300))

    outname = outid + '_ccs_mean_err.pdf'

    plt.savefig(os.path.join(dirpath, outname))
    plt.close()


def plot_ccs_distrbution_clust(dataframe, dirpath, outid):
    """
    plot ccs distribution according to cluster
    :param dataframe: dataframe
    :param dirpath: directory to save output
    :return: void.
    """

    outname = outid + '_ccs_distribution.pdf'

    with PdfPages(os.path.join(dirpath, outname)) as pdf:

        ccs_all = dataframe['CCS_TM'].values

        plt.hist(ccs_all, bins='auto', color='dodgerblue')
        plt.xlabel('IMPACT_CCS_TM')
        plt.ylabel('Count')
        plt.title('ALL')
        pdf.savefig()
        plt.close()

        unique_clust_num = np.unique(dataframe['clust_num'])
        for ind, clust_num in enumerate(unique_clust_num):

            df_clust = dataframe[dataframe['clust_num'] == clust_num]
            clust_ccs_tm = df_clust['CCS_TM'].values

            plt.hist(clust_ccs_tm, bins='auto', color='dodgerblue')
            plt.xlabel('IMPACT_CCS_TM')
            plt.ylabel('Count')
            plt.title('Clust_'+str(clust_num))
            pdf.savefig()
            plt.close()



def clust_pdb_ccs_stats(dataframe, dirpath, outid, linreg_attr=None):
    """
    write the stats from clust pdb ccs df
    :param dataframe: clust pdb df with ccs appended
    :param dirpath: directory to save output file
    :return: void
    """
    unique_clust_num = np.unique(dataframe['clust_num'])

    num_structs = []
    percent_structs = []

    mean_ccs_tm = []
    median_ccs_tm = []
    std_ccs_tm = []
    total_err_tm = []
    mean_ccs_pa = []
    median_ccs_pa = []
    std_ccs_pa = []


    for ind, clust_num in enumerate(unique_clust_num):

        total_num_structs = len(dataframe['clust_num'].values)
        df_clust = dataframe[dataframe['clust_num'] == clust_num]
        clust_ccs_pa = df_clust['CCS_PA'].values

        if linreg_attr:
            clust_ccs_tm = df_clust['CCS_TM'].values * linreg_attr[0] + linreg_attr[1]
        else:
            clust_ccs_tm = df_clust['CCS_TM'].values





        mean_ccs_tm.append(np.mean(clust_ccs_tm))
        median_ccs_tm.append(np.median(clust_ccs_tm))
        std_ccs_tm.append(np.std(clust_ccs_tm))


        if linreg_attr:
            cal_err = linreg_attr[2]*np.mean(clust_ccs_tm)/100
            total_err = sqrt(cal_err**2 + np.std(clust_ccs_tm)**2)
            total_err_tm.append(total_err)
        else:
            total_err_tm.append(np.std(clust_ccs_tm))


        mean_ccs_pa.append(np.mean(clust_ccs_pa))
        median_ccs_pa.append(np.median(clust_ccs_pa))
        std_ccs_pa.append(np.std(clust_ccs_pa))

        len_clust = len(df_clust['clust_num'].values)
        perc_struct = len_clust*100/total_num_structs

        num_structs.append(len_clust)
        percent_structs.append(perc_struct)




    output_string = ''

    if linreg_attr:
        linreg_attr_line = '#Slope,'+str(linreg_attr[0])+'\n#Intercept,'+str(linreg_attr[1])+'\n#RMSE(%),'+str(linreg_attr[2])+'\n'
        output_string += linreg_attr_line
    else:
        output_string = '#Slope,'+str(np.nan)+'\n#Intercept,'+str(np.nan)+'\n#RMSE(%),'+str(np.nan)+'\n'

    header = 'clust_num,num_structs,percent_structs,mean_ccs_tm,median_ccs_tm,std_ccs_tm,rmse_err,total_err,mean_ccs_pa,median_ccs_pa,std_ccs_pa\n'
    output_string += header

    if linreg_attr:
        rmse = linreg_attr[2]
    else:
        rmse = np.nan

    for num in range(len(unique_clust_num)):

        line = '{},{},{},{},{},{},{},{},{},{},{}\n'.format(unique_clust_num[num],
                                                     num_structs[num],
                                                     percent_structs[num],
                                                     mean_ccs_tm[num],
                                                     median_ccs_tm[num],
                                                     std_ccs_tm[num],
                                                     rmse,
                                                     total_err_tm[num],
                                                     mean_ccs_pa[num],
                                                     median_ccs_pa[num],
                                                     std_ccs_pa[num])
        output_string += line

    out_fname = outid + '_pdb_ccs_stats.csv'

    with open(os.path.join(dirpath, out_fname), 'w') as outfile:
        outfile.write(output_string)
        outfile.close()

    plot_ccs_mean_std(unique_clust_num, mean_ccs_tm, total_err_tm, dirpath, outid)



def gen_cluster_pdb_ccs_df(dirpath, cluster_pdb_file, impact_ccs_file):
    """
    generate cluster pdb ccs df
    :param dirpath: directory path
    :param cluster_pdb_file: cluster pdb file csv
    :param impact_ccs_file: impact ccs file csv
    :return: dataframe
    """

    clust_pdb_df = pd.read_csv(os.path.join(dirpath, cluster_pdb_file))
    impact_ccs_df = pd.read_csv(os.path.join(dirpath, impact_ccs_file))

    clust_pdb_ccs_pa = []
    clust_pdb_ccs_tm = []

    for index, clust_pdb_fpath in enumerate(clust_pdb_df['pdb_fname'].values):
        clust_pdb_fname = os.path.split(clust_pdb_fpath)[1]
        for ind, (pdb_name, ccs_pa, ccs_tm) in enumerate(zip(impact_ccs_df['filename'],
                                                             impact_ccs_df['CCS_PA'],
                                                             impact_ccs_df['CCS_TM'])):
            if clust_pdb_fname == pdb_name:
                clust_pdb_ccs_pa.append(ccs_pa)
                clust_pdb_ccs_tm.append(ccs_tm)


    clust_pdb_df['CCS_PA'] = clust_pdb_ccs_pa
    clust_pdb_df['CCS_TM'] = clust_pdb_ccs_tm

    return clust_pdb_df



def cluster_pdb_ccs_output(dirpath, cluster_pdb_file, impact_ccs_file, outid, linreg_attr=None):
    """
    write all the outputs
    :param dirpath: directory
    :param cluster_pdb_file: input cluster pdb ccs file
    :param impact_ccs_file: input impact ccs file
    :param outid: outid
    :return: void.
    """

    # write clust pdb ccs file
    outfname = outid + '_pdb_ccs.csv'
    clust_pdb_ccs_df = gen_cluster_pdb_ccs_df(dirpath, cluster_pdb_file, impact_ccs_file)
    clust_pdb_ccs_df.to_csv(os.path.join(dirpath, outfname), index=False)

    clust_pdb_ccs_stats(clust_pdb_ccs_df, dirpath, outid, linreg_attr=linreg_attr)

    plot_ccs_distrbution_clust(clust_pdb_ccs_df, dirpath, outid)

if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_nomod\rmsdmat"
    impact_ccs_file = 'repid_0_ccs.csv'
    clust_pdb_file = 'cluster_manual_pdb_files.csv'

    linreg_attr_mod = [0.739234, 495.4191,6.78335]
    linreg_attr_nomod = [0.919985, 148.8482, 6.963937]

    cluster_pdb_ccs_output(dirpath, clust_pdb_file, impact_ccs_file, outid='cluster_manual', linreg_attr=linreg_attr_nomod)