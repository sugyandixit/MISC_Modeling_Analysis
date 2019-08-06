# f_test_stats.py
import numpy as np
from scipy.stats import f_oneway, f

def between_group_variability(rmsd_clust):
    # rmsd_clust = np.array(rmsd_clust)
    group_means = np.array([np.mean(x) for x in rmsd_clust])
    overall_mean = np.mean(np.concatenate(rmsd_clust))
    len_groups = np.array([len(x) for x in rmsd_clust])
    expl_var_ = 0
    df_n = len(group_means) - 1
    for ind, (mean_group, num_group) in enumerate(zip(group_means, len_groups)):
        expl_var = (num_group * (mean_group - overall_mean)**2)/df_n
        expl_var_ += expl_var
    return expl_var_, df_n

def within_group_variability(rmsd_clust):
    group_means = np.array([np.mean(x) for x in rmsd_clust])
    overall_sample_size = len(np.concatenate(rmsd_clust))
    var_sum_ = 0
    df_d = overall_sample_size - len(group_means)
    for ind, (mean_group, rmsd_group) in enumerate(zip(group_means, rmsd_clust)):
        var_sum = np.sum(np.square(np.subtract(rmsd_group, mean_group)))
        var_sum_ += var_sum
    unexpl_var_ = var_sum_/df_d
    return unexpl_var_, df_d

def f_test(rmsd_clust):
    expl_var, df_n = between_group_variability(rmsd_clust)
    unexpl_var, df_d = within_group_variability(rmsd_clust)
    f_val = expl_var/unexpl_var
    pdf_ = f.pdf(f_val, df_n, df_d)
    # cdf_ = f.cdf(f_val, df_n, df_d)
    return expl_var, unexpl_var, f_val, pdf_