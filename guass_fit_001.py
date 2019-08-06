import os
import  numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gauss_func(x, y0, amp, xc, sigma):
    rxc = ((x - xc) ** 2) / (2 * (sigma ** 2))
    y = y0 + amp * (np.exp(-rxc))
    return y

def get_files(dirpath, end_file_id):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith(end_file_id)]
    return files

def histdata(data):
    hist, binedges = np.histogram(data, bins='auto', density=True)
    bincenters = (binedges[:-1]+binedges[1:])/2
    return hist, bincenters

def estimate_params_gauss(xdata, ydata, y0=0):
    center = np.median(xdata)
    amp = np.max(ydata)
    width = 0.05*len(xdata)
    return y0, amp, center, width

if __name__=='__main__':
    dirpath = r"C:\Users\sugyan\Documents\Papers\021317_UtilizingPeakwidth\_Submission_to_JASMS\ConstTempEnergyData"
    end_file_id = 'energy_001.csv'
    outline = ''
    header ='num,fname,y0,amp,center,sigma\n'
    outline += header
    files = get_files(dirpath, end_file_id)
    for ind, file in enumerate(files):
        fname = str(file)
        df = pd.read_csv(os.path.join(dirpath, file), sep='\s+', header=None)
        energy = df.iloc[100:,1]
        ydata, xdata = histdata(energy)
        std_energy = np.std(energy)
        y0, amp, center, width = estimate_params_gauss(xdata, ydata)
        popt, pcov = curve_fit(gauss_func, xdata, ydata, p0=(y0, amp, center, width))
        gauss_fit_pars = ','.join([str(x) for x in popt])
        line = '{},{},{}\n'.format(ind+1, fname, gauss_fit_pars)
        outline += line
        y_fit = gauss_func(xdata, *popt)
        plt.scatter(xdata, ydata, marker='o', color='black')
        plt.hist(energy, bins='auto', density=True, color='black', alpha=0.5)
        plt.plot(xdata, y_fit, linestyle='--', marker=None, color='red')
        plt.savefig(os.path.join(dirpath, fname+'_gaussfit.png'))
        plt.close()
        print('heho')

    outfile = end_file_id.split('.')[0] + '_gaussfit.csv'
    with open(os.path.join(dirpath, outfile), 'w') as gaussout:
        gaussout.write(outline)
        gaussout.close()