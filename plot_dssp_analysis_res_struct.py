import os
import pandas as pd
import matplotlib.pyplot as plt


def get_files(dirpath, endid):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith(endid)]
    return files

def plot_sec_struct_res(df, dirpath, outid, image_type='.pdf'):


    #plot_helix_beta_loop

    # plt.plot(df['res_num'].values, df['frac_helix'].values, ls='-', marker='o', color='magenta', label='helix')
    plt.plot(df['res_num'].values, df['frac_310_helix'].values, ls='-', marker='o', color='green', label='310_helix')
    plt.plot(df['res_num'].values, df['frac_alpha_helix'].values, ls='-', marker='o', color='dodgerblue', label='alpha_helix')
    plt.plot(df['res_num'].values, df['frac_pi_helix'].values, ls='-', marker='o', color='red', label='pi_helix')
    plt.plot(df['res_num'].values, df['frac_beta_strand'].values, ls='-', marker='o', color='cyan', label='beta_strand')
    plt.plot(df['res_num'].values, df['frac_loop'].values, ls='-', marker='o', color='wheat', label='loop')
    plt.xlabel('RES')
    plt.ylabel('Fraction')
    plt.legend(loc='best', fontsize='small')


    outname = outid + image_type

    plt.savefig(os.path.join(dirpath, outname))
    plt.close()


if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_nomod\rmsdmat"
    dssp_analysis_res_files = get_files(dirpath, endid='dssp_analysis_res.csv')
    for file in dssp_analysis_res_files:
        outid = str(file).split('.csv')[0]
        df = pd.read_csv(os.path.join(dirpath, file))
        plot_sec_struct_res(df, dirpath, outid, image_type='.png')