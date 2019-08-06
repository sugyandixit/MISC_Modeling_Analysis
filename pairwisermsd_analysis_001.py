import os
import numpy as np
import matplotlib.pyplot as plt

def get_rmsd_mat_file(dirpath):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('.npy')]
    return files

def read_rmsd_data(dirpath, file):
    rmsddata = np.load(os.path.join(dirpath, file))
    return rmsddata

def write_stats(rmsddata, dirpath, temp, fname_id):
    mean = np.mean(rmsddata)
    median = np.median(rmsddata)
    std = np.std(rmsddata)
    fname = fname_id+'_'+temp+'_statout.csv'
    with open(os.path.join(dirpath, fname), 'w') as statout:
        header = '#temp, mean, median, std\n'
        statout.write(header)
        line = '{}, {}, {}, {}\n'.format(temp, mean, median, std)
        statout.write(line)
    statout.close()

def plot_rmsd_hist(rmsddata, dirpath, temp, fname_id):
    fname = fname_id + '_'+temp+'_hisplot.pdf'
    plt.hist(rmsddata, bins='auto', density=True)
    plt.xlabel('Pairwise RMSD (Angstrom)')
    plt.ylabel('Count')
    plt.savefig(os.path.join(dirpath, fname))
    plt.close()


if __name__=='__main__':
    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\AcAlaPeptide_ClusterAnalysis\helix_018"
    dir_p, dir_n = os.path.split(dirpath)
    rmsdfiles = get_rmsd_mat_file(dirpath)
    for file in rmsdfiles:
        rmsddata = read_rmsd_data(dirpath, file)
        file_chars_dot = str(file).split('.')[0]
        temp_ids = file_chars_dot.split('_structs_')[-1]
        fnameid = dir_n+'_rmsdmatrix'
        write_stats(rmsddata, dirpath, temp_ids, fnameid)
        plot_rmsd_hist(rmsddata, dirpath, temp_ids, fnameid)
