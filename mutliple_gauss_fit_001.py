import os
import numpy as np
from scipy.stats import linregress
from scipy.optimize import leastsq
from math import sqrt, pi
import matplotlib.pyplot as plt

def gauss_func(x, params):
    y0, amp, xc, sigma = params
    rxc = ((x - xc) ** 2) / (2 * (sigma ** 2))
    y = y0 + amp * (np.exp(-rxc))
    return y

def multiple_gauss_func(num, x, params):
    params_split = np.split(params, num)
    y_store = []
    for ind in range(num):
        y = gauss_func(x, params_split[ind])
        y_store.append(y)
    y_fit = np.sum(y_store, axis=0)
    y_store = np.asarray(y_store)
    return y_fit, y_store

def histdata(data):
    hist, binedges = np.histogram(data, bins='auto', density=True)
    bincenters = (binedges[:-1]+binedges[1:])/2
    return hist, bincenters

def min_function_gausfit(params, num, y, x):
    yfit = multiple_gauss_func(num, x, params)[0]
    cost = np.square(np.subtract(y, yfit))
    return cost

def residual_func(ydata, yfit):
    res = ydata - yfit
    return res

def variance_from_residual(residual):
    y = np.sum(np.square(residual))/(len(residual)-2)
    return y

def write_gausfits(xdata, ydata, yfit, yfit_components, residual, data_fname, dirpath):
    num_components = np.arange(1, len(yfit_components)+1)
    ycomp_label = ','.join('y_'+str(x) for x in num_components)
    header = '# xdata, ydata, ' + ycomp_label + ', y_cumulative, residual'
    output_string = ''
    output_string += header +'\n'
    for index in range(len(yfit_components[0])):
        x_data = str(xdata[index])
        y_data = str(ydata[index])
        joined_y_components = ','.join([str(x) for x in yfit_components[:, index]])
        y_fit_cumulative = str(yfit[index])
        residual_ = str(residual[index])
        line1 = '{}, {}, {}, {}, {}\n'.format(x_data, y_data, joined_y_components, y_fit_cumulative, residual_)
        output_string += line1

    with open(os.path.join(dirpath, data_fname+'_gaussfit_output.csv'), 'w') as gausfitout:
        gausfitout.write(output_string)

def gauss_weights(num, gauss_params):
    gauss_params = np.split(gauss_params, num)
    peak_area = []
    for ind, array in enumerate(gauss_params):
        gauss_pk_area = array[1] * array[3] * sqrt(2*pi)
        peak_area.append(gauss_pk_area)
    peak_area = np.asarray(peak_area)
    gauss_wts = []
    for ind in range(len(peak_area)):
        wt = peak_area[ind]/np.sum(peak_area)
        gauss_wts.append(wt)
    gauss_wts = np.asarray(gauss_wts)
    return peak_area, gauss_wts

def wt_mean_var_mixture_stat(num, gauss_params):

    gauss_wts = gauss_weights(num, gauss_params)[1]
    gauss_params = np.split(gauss_params, num)
    wt_mu = 0
    wt_var = 0
    for index, array in enumerate(gauss_params):
        mu = array[2]
        # var = array[3]
        wt = gauss_wts[index]
        wt_mu0 = wt*mu
        wt_mu = wt_mu + wt_mu0
    for ind, array in enumerate(gauss_params):
        mu = array[2]
        var = array[3]
        wt = gauss_wts[ind]
        wt_var0 = wt*(var+np.square(mu - wt_mu))
        wt_var = wt_var + wt_var0
    return wt_mu, wt_var

def write_gausfit_statout(num, gauss_params, gauss_area, gauss_wt, wt_mean, wt_var, variance, lnregression, data_fname, dirpath):
    # num_components = np.arange(1, len(gauss_params)+1)
    header1 = '# component, y0, amp, xc, w, area, weights\n'
    output_string = ''
    output_string += header1
    gauss_params = np.split(gauss_params, num)
    for index, array in enumerate(gauss_params):
        join_gauss_params = ','.join([str(x) for x in array])
        line1 = '{}, {}, {}, {}\n'.format(str(index), join_gauss_params, gauss_area[index], gauss_wt[index])
        output_string += line1
    header2 = '# combined_stats, wt_mean, wt_var, slope, variance, r2\n'
    output_string += header2
    line2 = '{}, {}, {}, {}, {}, {}\n'.format('', wt_mean, wt_var, lnregression[0], variance, lnregression[2]**2)
    output_string += line2

    with open(os.path.join(dirpath, data_fname + '_gausfit_statout.csv'), 'w') as statout:
        statout.write(output_string)



def plot_residual(x, residual, data_fname, dirpath):
    plt.scatter(x, residual, color='black')
    plt.xlabel('Pairwise rmsddistribution')
    plt.ylabel('residual')
    plt.savefig(os.path.join(dirpath, data_fname+'_gaussfit_residual.pdf'))
    plt.close()

def plot_gauss_fits(xdata, rmsddata, yfit, yfit_component, data_fname, dirpath):
    colors = ['red', 'blue', 'magenta', 'green', 'orange', 'lightblue']
    plt.hist(rmsddata, bins='auto', density=True, color='black', alpha=0.5)
    plt.plot(xdata, yfit, color='black', ls='--')
    for ind, arr in enumerate(yfit_component):
        plt.plot(xdata, arr, color=colors[ind])
    plt.savefig(os.path.join(dirpath, data_fname+'_gaussfit.pdf'))
    plt.close()

if __name__=='__main__':
    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\AcAlaPeptide_ClusterAnalysis\helix_019"
    dirp, dirn = os.path.split(dirpath)
    data_fname = "rmsdonlymat_1000_structs_300.npy"
    output_fname = dirn + '_' + data_fname.split('.')[0]
    rmsd_data = np.load(os.path.join(dirpath, data_fname))
    ydata, xdata = histdata(rmsd_data)

    num_gauss = 3
    # params_1 = [0, 0.44, 1.4, 0.3]
    # params_2 = [0, 0.39, 2.4, 0.5]
    # params_3 = [0, 0.27, 3.8, 0.2]
    # params_4 = [0, 0.1, 4.5, 0.2]

    params_1 = [0, 0.9, 1.1, 0.26]
    params_2 = [0, 0.9, 1.2, 0.4]
    params_3 = [0, 0.3, 1.8, 0.3]
    params_4 = [0, 0.05, 3.0, 0.3]
    params_4_ = [0, 0.08, 3.0, 0.5]
    params_5 = [0., 0.005, 4.8, 0.2]
    params_6 = [0., 0.07, 4.7, 0.5]

    params_list = [params_1, params_2, params_3, params_4, params_5, params_6]
    params = []
    for pars in range(num_gauss):
        params.append(params_list[pars])
    params = np.concatenate(params)

    sol, flag = leastsq(min_function_gausfit, params, args=(num_gauss, ydata, xdata), maxfev=10000, ftol=0.0001)

    sol_ysub = []
    for arr in np.split(sol, num_gauss):
        arr_ = [0, arr[1], arr[2], arr[3]]
        sol_ysub.append(arr_)

    sol_ysub = np.asarray(sol_ysub)
    sol_ysub = np.concatenate(sol_ysub)

    # sol_ysub = sol



    yfit, yfit_component = multiple_gauss_func(num_gauss, xdata, sol_ysub)

    area, gauss_wts = gauss_weights(num_gauss, sol_ysub)
    wt_mean, wt_var = wt_mean_var_mixture_stat(num_gauss, sol_ysub)

    lnreg = linregress(ydata, yfit)
    residual = residual_func(ydata, yfit)
    variance = variance_from_residual(residual)

    write_gausfits(xdata, ydata, yfit, yfit_component, residual, output_fname, dirpath)
    write_gausfit_statout(num_gauss, sol_ysub, area, gauss_wts, wt_mean, wt_var, variance, lnreg, output_fname, dirpath)

    plot_gauss_fits(xdata, rmsd_data, yfit, yfit_component, output_fname, dirpath)
    plot_residual(xdata, residual, output_fname, dirpath)
