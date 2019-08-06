import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as hclust
import matplotlib.pyplot as plt


fname = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\jac_dist_ct_map_10.csv"
data = np.genfromtxt(fname, delimiter=',', skip_header=1, dtype='str')
jac_dist_numpy = np.array(data[:, 2], dtype='float')
jac_dist = pd.read_csv(fname, sep=',', usecols=['jac_dist'])
jac_dist = jac_dist.jac_dist.values
# jac_dist = np.concatenate(jac_dist.values)
# jac_dist = jac_dist.values.reshape(-1, 1)

jac_dist_sqform = squareform(jac_dist)
jac_dist_triu = np.triu(jac_dist_sqform)
plt.imshow(jac_dist_triu, cmap='binary', interpolation='nearest', origin='lower')
plt.show()
plt.close()

linkage = hclust.linkage(jac_dist_sqform, method='average')
plt.figure()
dn = hclust.dendrogram(linkage, distance_sort='ascending')
plt.show()
plt.close()