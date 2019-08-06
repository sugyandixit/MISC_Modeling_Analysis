# creates charge permutation in a sequence
# author: Suggie
# date: 02/04/19

import os
import itertools
import numpy as np
import random
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy.stats import zscore
from matplotlib_venn import venn2, venn3

def read_sequence(seq_file):
    """
    read the sequence file and append sequence into an array
    :param seq_file: sequence file
    :return: sequence array (numpy)
    """
    seq = []
    sequ_f = open(seq_file, 'r')
    sequ_f = sequ_f.read().splitlines()
    for line in sequ_f:
        chars = line.split()
        for res in chars:
            seq.append(res.upper())
    seq = np.array(seq)
    return seq

def divide_sequence(seq, num=2):
    """
    divide the sequence to equal or near equal sub arrays indicated by the number
    :param seq: sequence list
    :param num: number of divisions
    :return: sub array
    """
    split_seq_array = np.array_split(seq, num)
    return split_seq_array


def seq_encoder(sequence, nter=True, basic_res=['LYS', 'ARG'], acidic_res=['ASP', 'GLU']):
    """
    anotate residue by integers based upon whether they're neutral (0) / basic(1) / acidic(-1)
    :param sequence: sequence list
    :param basic_res: res that are basic
    :param acidic_res: res that are acidic
    :return: return sequence as ints
    """
    seq_int = []
    for ind, res in enumerate(sequence):
        if res == basic_res[0] or res == basic_res[1]:
            seq_int.append(1)
        elif res == acidic_res[0] or res == acidic_res[1]:
            seq_int.append(-1)
        else:
            seq_int.append(0)
    if nter == True:
        seq_int = [1] + seq_int
    else:
        seq_int = seq_int
    seq_int = np.array(seq_int)
    return seq_int

def gen_charge_dict_from_seq(seq_ints):
    charge = np.sum(seq_ints)
    basic_label = 1
    acidic_label = -1
    neutral_label = 0
    basic_sites_inds = []
    acidic_sites_inds = []
    neutral_sites_inds = []
    basic_num = 0.0
    acidic_num = 0.0
    neutral_num = 0.0
    for ind, seq in enumerate(seq_ints):
        if seq == basic_label:
            basic_num += 1
            basic_sites_inds.append(ind)
        if seq == acidic_label:
            acidic_num += 1
            acidic_sites_inds.append(ind)
        if seq == neutral_label:
            neutral_num += 1
            neutral_sites_inds.append(ind)
    charge_dict = dict()
    # charge_dict['seq_int'] = seq_ints
    charge_dict['num_charge'] = int(charge)
    charge_dict['num_basic_sites'] = int(basic_num)
    charge_dict['num_acidic_sites'] = int(acidic_num)
    charge_dict['num_neutral_sites'] = int(neutral_num)
    charge_dict['basic_sites_indices'] = basic_sites_inds
    charge_dict['acidic_sites_indices'] = acidic_sites_inds
    charge_dict['neutral_sites_indices'] = neutral_sites_inds
    return charge_dict


def generate_charge_site_combinations(sample, size, random_sample_size=50, max_sample_size=1000, filter_distance=None):
    comb_list = []
    for index in range(max_sample_size):
        if index < random_sample_size:
            comb = sorted(random.sample(sample, size))
            if None == filter_distance:
                comb_list.append(np.array(comb))
            else:
                charge_site = filter_by_charge_dist(np.array(comb), filter_distance)
                if charge_site is not None:
                    comb_list.append(charge_site)
    return comb_list


def filter_by_charge_dist(charge_sites, dist=3):
    """
    filter the charge sites list by discarding any that have charges within a distance of 3
    :param charge_sites: charge sites array
    :param dist: distance at which the charge sites can be located
    :return: list of charge sites array
    """
    # out_charge_site = []
    if len(charge_sites) == 1:
        out_charge_site = charge_sites
    else:
        for ind in range(len(charge_sites)):
            if (ind+1) < len(charge_sites):
                dist_res = charge_sites[ind+1] - charge_sites[ind]
                if dist_res <= dist:
                    return None
                else:
                    out_charge_site = charge_sites
    return out_charge_site



def charge_permutation(sequence, nter=False, charge_assign_num=4, random_sample_size=10):
    """
    perform charge permutation with either basic sites or acidic sites
    :param sequence: sequence list
    :param charge_num: number of basic sites to be indicated as having charge
    :return: list of permutation of sequence with basic sites indicated as having charge
    """

    charge_perm_obj = ChargePermutation()
    charge_perm_obj.sequence = sequence
    charge_perm_obj.nter = nter
    charge_perm_obj.assign_charge = charge_assign_num

    seq_label = seq_encoder(sequence, nter=nter)
    seq_charge_dict = gen_charge_dict_from_seq(seq_label)

    charge_perm_dict_list =[]
    for num_negative_charge in range(1, seq_charge_dict['num_acidic_sites']+1):
        for num_positive_charge in range(1, seq_charge_dict['num_basic_sites']+1):
            enum_charge = num_positive_charge - num_negative_charge
            if enum_charge == charge_assign_num:
                charge_perm_dict = dict()
                charge_perm_dict['num_pos_charge'] = num_positive_charge
                charge_perm_dict['num_neg_charge'] = num_negative_charge
                charge_perm_dict_list.append(charge_perm_dict)

    for charge_perms in charge_perm_dict_list:
        # print('pos')
        pos_perm_indices = generate_charge_site_combinations(seq_charge_dict['basic_sites_indices'], charge_perms['num_pos_charge'], filter_distance=3, random_sample_size=random_sample_size)
        # print('neg')
        neg_perm_indices = generate_charge_site_combinations(seq_charge_dict['acidic_sites_indices'], charge_perms['num_neg_charge'], filter_distance=3, random_sample_size=random_sample_size)
        for neg_perm_ind in neg_perm_indices:
            for pos_perm_ind in pos_perm_indices:
                temp_seq_int = np.zeros_like(seq_label)
                for neg_res_ind in neg_perm_ind:
                    temp_seq_int[neg_res_ind] = -1
                for pos_res_ind in pos_perm_ind:
                    temp_seq_int[pos_res_ind] = 1
                temp_seq_charge_dict = gen_charge_dict_from_seq(temp_seq_int)
                charge_perm_obj.seq_ints.append(temp_seq_int)
                charge_perm_obj.charge_dict.append(temp_seq_charge_dict)

    return charge_perm_obj



class ChargePermutation(object):
    """
    stores the charge permutation information. Includes the num_positive and num_negative charge, pos_charge_index,
    neg_charge_index, sequence in charge, sequence information, Nter=true or False
    """

    def __init__(self):
        """
        initialize the object. Requires the sequence, nter bool, and charges to be assigned.
        """
        self.sequence = None
        self.nter = None
        self.assign_charge = None
        self.charge_dict = []
        self.seq_ints = []



def plot_charge_dist(list_seq_ints, dirpath, fname):
    cmap = ListedColormap(['red', 'lightgrey', 'blue'])
    bounds = [-2,0,1,2]
    norm = BoundaryNorm(bounds, cmap.N)
    seq_charge_stack = np.vstack(list_seq_ints)
    plt.pcolor(seq_charge_stack, cmap=cmap, norm=norm)
    plt.savefig(os.path.join(dirpath, fname+'.png'))
    plt.close()
    print('heho')



def charge_comb_with_sequ_div(sequence, div_num=4, nter_list=[True, False, False, False], charge_per_seg=3, random_sample_size_per_seg=2):
    seq_div = divide_sequence(sequence, num=div_num)
    charge_comb_seq_int_list = [[] for _ in range(div_num)]
    for ind, seq in enumerate(seq_div):
        charge_comb = charge_permutation(seq, nter=nter_list[ind], charge_assign_num=charge_per_seg, random_sample_size=random_sample_size_per_seg)
        seq_int = charge_comb.seq_ints
        charge_comb_seq_int_list[ind].append(seq_int)

    product_list = []
    for prod in itertools.product(*charge_comb_seq_int_list[0], *charge_comb_seq_int_list[1], *charge_comb_seq_int_list[2], *charge_comb_seq_int_list[3]):
        seq_int = np.concatenate(prod)
        product_list.append(seq_int)

    return product_list


def get_basic_acidic_charge_ind(seq_int_list):
    charge_sum_array = np.sum(seq_int_list, axis=0)
    basic_ind_arr = []
    basic_charge_arr = []
    acidic_ind_arr = []
    acidic_charge_arr = []
    for num, charge in enumerate(charge_sum_array):
        if charge > 0:
            basic_charge_arr.append(charge)
            basic_ind_arr.append(num)
        if charge < 0:
            acidic_charge_arr.append(charge)
            acidic_ind_arr.append(num)

    return basic_charge_arr, basic_ind_arr, acidic_charge_arr, acidic_ind_arr


def plot_z_score(arr, arr_label, dirpath, fname):
    z_score = zscore(arr)
    y_pos = np.arange(len(z_score))
    plt.bar(y_pos, z_score, align='center', alpha=0.7, color='black')
    plt.xticks(y_pos, arr_label)
    plt.xlabel('res_index')
    plt.ylabel('Zscore')
    plt.savefig(os.path.join(dirpath, fname+'_zscore.png'))
    plt.close()


def plot_2venn_diagram(arr_list, arr_label_list, dirpath, fname):
    v = venn2(arr_list, set_labels=arr_label_list)
    plt.savefig(os.path.join(dirpath, fname+'_venndiagram.png'))
    plt.close()

def plot_3venn_diagram(arr_list, arr_label_list, dirpath, fname):
    v = venn3(arr_list, set_labels=arr_label_list)
    plt.savefig(os.path.join(dirpath, fname+'_venndiagram.png'))
    plt.close()


if __name__=='__main__':
    seq_file = r"C:\Users\sugyan\Documents\MembraneMDfiles\ChargePermutation\HSA_D12_3letter.txt"
    dpath = os.path.split(seq_file)[0]
    sequence = read_sequence(seq_file)
    seq_int = seq_encoder(sequence, nter=True)
    charge_dict_sequence = gen_charge_dict_from_seq(seq_int)
    charge_comb = charge_permutation(sequence, nter=True, charge_assign_num=12, random_sample_size=40)
    full_charge_comb_summary = get_basic_acidic_charge_ind(charge_comb.seq_ints)
    # plot_charge_dist(charge_comb.seq_ints, dirpath=dpath, fname='full_seq_charge_placement')
    div_seq_int_list = charge_comb_with_sequ_div(sequence, div_num=4, nter_list=[True, False, False, False], charge_per_seg=3, random_sample_size_per_seg=2)
    # plot_charge_dist(div_seq_int_list, dirpath=dpath, fname='div_seq_charge_placement')
    div_charge_comb_summary = get_basic_acidic_charge_ind(div_seq_int_list)

    plot_z_score(full_charge_comb_summary[0], full_charge_comb_summary[1], dirpath=dpath, fname='full_seq_basic_sites')
    plot_z_score(full_charge_comb_summary[2], full_charge_comb_summary[3], dirpath=dpath, fname='full_seq_acidic_sites')

    plot_z_score(div_charge_comb_summary[0], div_charge_comb_summary[1], dirpath=dpath, fname='div_seq_basic_sites')
    plot_z_score(div_charge_comb_summary[2], div_charge_comb_summary[3], dirpath=dpath, fname='div_seq_acidic_sites')


    basic_venn_list_full = [set(charge_dict_sequence['basic_sites_indices']), set(full_charge_comb_summary[1])]
    acidic_venn_list_full = [set(charge_dict_sequence['acidic_sites_indices']), set(full_charge_comb_summary[3])]
    venn_label_list_full = ['All', 'full_seq_charge_comb']

    basic_venn_list_div = [set(charge_dict_sequence['basic_sites_indices']), set(div_charge_comb_summary[1])]
    acidic_venn_list_div = [set(charge_dict_sequence['acidic_sites_indices']), set(div_charge_comb_summary[3])]
    venn_label_list_div = ['All', 'div_seq_charge_comb']

    plot_2venn_diagram(basic_venn_list_full, venn_label_list_full, dirpath=dpath, fname='full_basic_sites')
    plot_2venn_diagram(acidic_venn_list_full, venn_label_list_full, dirpath=dpath, fname='full_acidic_sites')

    plot_2venn_diagram(basic_venn_list_div, venn_label_list_div, dirpath=dpath, fname='div_basic_sites')
    plot_2venn_diagram(acidic_venn_list_div, venn_label_list_div, dirpath=dpath, fname='div_acidic_sites')

    basic_venn_list = [set(charge_dict_sequence['basic_sites_indices']), set(full_charge_comb_summary[1]), set(div_charge_comb_summary[1])]
    acidic_venn_list = [set(charge_dict_sequence['acidic_sites_indices']), set(full_charge_comb_summary[3]), set(div_charge_comb_summary[3])]
    venn_label_list = ['All', 'full_seq_charge_comb', 'div_seq_charge_comb']

    plot_3venn_diagram(basic_venn_list, venn_label_list, dirpath=dpath, fname='basic_sites')
    plot_3venn_diagram(acidic_venn_list, venn_label_list, dirpath=dpath, fname='acidic_sites')










    # dividing the sequence
    # seq_div = divide_sequence(sequence, num=4)
    # charge_comb_obj_list = []
    # nter_list = [True, False, False, False]
    # for ind, seq in enumerate(seq_div):
    #     charge_comb = charge_permutation(seq, nter=nter_list[ind], charge_assign_num=3, random_sample_size=2)
    #     charge_comb_obj_list.append(charge_comb)
    print('heho')