# plot three dimensional landscape
# need similarity square matrix and scoring/energy matrix

import os
import numpy as np
from sklearn import manifold
from scipy.interpolate import griddata
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation



def plot_mds_space(pos, dirpath):
    plt.scatter(pos[1:, 0], pos[1:, 1], color='black', label='MDS')
    plt.savefig(os.path.join(dirpath, 'mds_metricTrue_.png'))
    plt.close()


def compute_mds(dist_mat_sqform, dirpath):
    seed = np.random.RandomState(seed=3)
    mds = manifold.MDS(n_components=2, max_iter=300, eps=1e-6, random_state=seed, dissimilarity='precomputed',
                       n_jobs=-1, n_init=4, verbose=2)
    pos = mds.fit_transform(dist_mat_sqform)
    plot_mds_space(pos, dirpath)
    return pos


def generate_grid_data(dist_mat_sqform, scoring_matrix, dirpath, num_points=100, savgol=True, window_savgol=25):
    pos = compute_mds(dist_mat_sqform, dirpath)

    x, y = pos[:, 0], pos[:, 1]
    x = x.T

    xi = np.linspace(min(x), max(x), num_points)
    yi = np.linspace(min(y), max(y), num_points)
    zi = griddata((x, y), scoring_matrix, (xi[None, :], yi[:, None]), method='nearest', fill_value=0)

    xim, yim = np.meshgrid(xi, yi)

    zim = []

    if savgol == True:
        norm = savgol_filter(zi, window_savgol, 2, axis=0)
        norm = savgol_filter(norm, window_savgol, 2, axis=1)
        zim = norm

    grid_data = [xim, yim, zim]

    return grid_data

def plot_landscape(dist_mat_sqform, scoring_matrix, dirpath, elev=13, azim=60, num_points=100, savgol=True, window_savgol=25):
    grid_data = generate_grid_data(dist_mat_sqform, scoring_matrix, dirpath, num_points=100, savgol=True,
                                   window_savgol=window_savgol)
    fig = plt.figure()
    ax = Axes3D(fig)

    m = cm.ScalarMappable(cmap='jet')
    m.set_array(grid_data[2])

    ax.plot_surface(grid_data[0], grid_data[1], grid_data[2], rstride=1, cstride=1, cmap = 'jet')
    ax.view_init(elev=elev, azim=azim)
    plt.savefig(os.path.join(dirpath, '3dsurface.png'))
    plt.close()

def animate_init_func(dist_mat_sqform, scoring_matrix, dirpath, elev=13, azim=140, num_points=100, savgol=True, window_savgol=25):
    grid_data = generate_grid_data(dist_mat_sqform, scoring_matrix, dirpath, num_points=100, savgol=True,
                                   window_savgol=25)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(grid_data[0], grid_data[1], grid_data[2], rstride=1, cstride=1, cmap='jet')
    return ()

def animate_func(i, elev=13, azim=140):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.view_init(elev=elev, azim=i+azim)
    return ()


def animation_landscape(dist_mat_sqform, scoring_matrix, dirpath, elev=13, azim=140, num_points=100, savgol=True, window_savgol=25):
    fig = plt.figure()
    ax = Axes3D(fig)
    init_func = animate_init_func(dist_mat_sqform, scoring_matrix, dirpath, elev=13, azim=140, num_points=100, savgol=True, window_savgol=25)

    anim = animation.FuncAnimation(fig, animate_func, init_func = animate_init_func, frames=720, interval=10, blit=True)
    anim.save('test_animation_001.mp4', fps=30, extra_args=['-vcodec', 'libx264'])