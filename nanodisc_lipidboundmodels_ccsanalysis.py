import os
import numpy as np
import matplotlib.pyplot as plt

from ccs_read import read_impact_ccs_file as read_impact
from ccs_read import conv_imp_to_csv

def plot_delta_ccs_vs_lipids(num_arr, ccs_pa_mean_arr, ccs_nolipid, dirpath):
    ccs_init = ccs_nolipid
    num_store = []
    delta_ccs = []
    for ind, num in enumerate(num_arr):
        d_ccs = ccs_pa_mean_arr[ind] - ccs_init
        delta_ccs.append(d_ccs)
        num_store.append(num)
        ccs_init = ccs_pa_mean_arr[ind]

    plt.scatter(num_store, delta_ccs, color='black')
    plt.xlabel('Num of lipids')
    plt.ylabel('Delta CCS')
    plt.savefig(os.path.join(dirpath, 'delta_ccs_vs_lipids.pdf'))
    plt.close()

def plot_ccs_vs_lipids(num_arr, ccs_pa_mean_arr, ccs_pa_std_arr, dirpath, ccs_nolipid=True):
    if ccs_nolipid:
        plt.scatter(0, ccs_nolipid, color='black')
    plt.errorbar(num_arr, ccs_pa_mean_arr, yerr=ccs_pa_std_arr, color='black', marker='o', ls='none')
    plt.xlabel('Num of lipids')
    plt.ylabel('CCS ($A^2$)')
    plt.savefig(os.path.join(dirpath, 'ccs_vs_lipids.pdf'))
    plt.close()

def plot_ccs_histogram(ccs, dirpath, fname):
    plt.hist(ccs, bins='auto', color='lightblue')
    plt.xlabel('CCS_IMP_TM')
    plt.ylabel('Frequency')
    plt.savefig(os.path.join(dirpath, fname+'_.pdf'))
    plt.close()

def save_ccs_numlipids_txt(num_lipids, ccs_tm_mean, ccs_tm_std, ccs_pa_mean, ccs_pa_std, dirpath):
    fname = 'numlipids_ccs_out.csv'
    with open(os.path.join(dirpath, fname), 'w') as saveccs:
        header = '{}, {}, {}, {}, {}\n'.format('num_lipids', 'ccs_tm_mean', 'ccs_tm_std', 'ccs_pa_mean', 'ccs_pa_std')
        saveccs.write(header)
        for ind, num in enumerate(num_lipids):
            line = '{}, {}, {}, {}, {}\n'.format(num, ccs_tm_mean[ind], ccs_tm_std[ind], ccs_pa_mean[ind], ccs_pa_std[ind])
            saveccs.write(line)
    saveccs.close()






if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\charmm_gui_iapp_nanodisc\iappdimer_peripheral_2\step_9_simanneal\numlipids_1\step9.4_low"
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('.ccsout')]

    num_lipid_all = []
    ccs_pa_all = []
    ccs_tm_all = []

    # num_lipid = []

    for file in files:
        fname, ccs_pa, ccs_tm = read_impact(file, dirpath)
        ccs_pa = np.asarray(ccs_pa, dtype='float')
        # print(ccs_pa)
        ccs_tm = np.asarray(ccs_tm, dtype='float')
        frame_num = []
        for item in fname:
            chars = item.split('_')
            frame_num.append(int(chars[-2]))
        conv_imp_to_csv(fname, ccs_pa, ccs_tm, dirpath, frame=frame_num)

    stop



    for file in files:
        fname2, ccs_pa, ccs_tm = read_impact(file, dirpath)
        ccs_pa_all.append(np.asarray(ccs_pa))
        ccs_tm_all.append(np.asarray(ccs_tm))
        fname = str(file).split('.')[0]
        numlipid = fname.split('_')[-1]
        num_lipid_all.append(numlipid)
        plot_ccs_histogram(np.asarray(ccs_tm), dirpath, fname)



    ccs_pa_mean_arr = []
    ccs_pa_std_arr = []
    ccs_tm_mean_arr = []
    ccs_tm_std_arr = []
    num_arr = []
    for ind, num_lipid_arr in enumerate(num_lipid_all):
        num = num_lipid_arr
        num_arr.append(num)
        ccs_pa_mean = np.mean(ccs_pa_all[ind])
        ccs_pa_std = np.std(ccs_pa_all[ind])
        ccs_pa_mean_arr.append(ccs_pa_mean)
        ccs_pa_std_arr.append(ccs_pa_std)
        ccs_tm_mean = np.mean(ccs_tm_all[ind])
        ccs_tm_std = np.std(ccs_tm_all[ind])
        ccs_tm_mean_arr.append(ccs_tm_mean)
        ccs_tm_std_arr.append(ccs_tm_std)

    sort_index = np.argsort(np.array(num_arr, dtype='int'))
    num_arr = np.array(num_arr)[sort_index]
    ccs_pa_mean_arr = np.array(ccs_pa_mean_arr)[sort_index]
    ccs_pa_std_arr = np.array(ccs_pa_std_arr)[sort_index]
    ccs_tm_mean_arr = np.array(ccs_tm_mean_arr)[sort_index]
    ccs_tm_std_arr = np.array(ccs_tm_std_arr)[sort_index]

    save_ccs_numlipids_txt(num_arr, ccs_tm_mean_arr, ccs_tm_std_arr, ccs_pa_mean_arr, ccs_pa_std_arr, dirpath)


    plot_ccs_vs_lipids(num_arr, ccs_tm_mean_arr, ccs_tm_std_arr, dirpath, ccs_nolipid=False)
    # plot_delta_ccs_vs_lipids(num_arr, ccs_tm_mean_arr, ccs_nolipid_tm, dirpath)

