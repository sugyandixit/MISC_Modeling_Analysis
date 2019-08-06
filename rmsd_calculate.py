# calculate rmsd

"""
Calculates rmsd between two pdb files

author:SDixit
date:03/17/17
"""

# import os
# import itertools
# import time
import numpy as np
# from scipy.spatial.distance import pdist, squareform

def get_coordinates(pdbfilepath):
    pdb = open(pdbfilepath, 'r+')
    pdb = pdb.read().splitlines()
    coor = []
    atoms = []
    for line in pdb:
        if line.startswith('ATOM'):
            x = line[30:38]
            y = line[38:46]
            z = line[46:54]
            atomid = line[13]
            coor.append(np.asarray([x, y, z], dtype='float'))
            atoms.append(np.asarray(atomid))
    coor = np.asarray(coor)
    atoms = np.asarray(atoms)
    assert(coor.shape[0] == atoms.size)
    return coor

def centroid(X):
    """
    Calculate the centroid from a vector set X
    :param X:
    :return:
    """
    centr = sum(X)/len(X)
    return centr

def rmsdcalc(P, Q):
    """
    Calculate the root mean square deviation from two sets of vectors V and W
    :param P:
    :param Q:
    :return:
    """
    D = len(P[0])
    N = len(P)
    rmsd = 0.0
    for p, q in zip(P, Q):
        rmsd += sum([(p[i] - q[i])**2.0 for i in range(D)])
    return np.sqrt(rmsd/N)



def kabsch(P, Q):
    """
    The optimal rotation matrix U is calculated and then used to rotate matrix P unto matrix Q so the minimum root mean
    square deviaion (RMSD) can be calculated

    Each set is represented as an N X D matrix, where D is the dimension in (x, y, and z (or more)) and N is number of
    points

    This algorithm works in three steps: a translation, the computation of a covariance matrix, and the computation of
    the optimal rotation matrix

    Translation - both sets of coordinates must be translated so that their centroids coincides with the origin of the
    coordinate system. Center of gravity to be (0, 0)

    Computation of the covariance matrix - calculate covariance matrix A with the two sets

    Computation of the optimal rotation matrix - first calculating the singular value decomposition (SVD) of the
    covariance matrix A. V, S, W = SVD(A). Check whether you need to correct the rotation matrix to ensure a right
    handed coordinate system. Then calculate optimal rotation matrix, U, as V. W

    :param P:
    :param Q:
    :return:
    """
    # create the centroid of P and Q which is the geometric center of a N-dimensional region and translate P and Q
    # onto that center
    Pc = centroid(P)
    Qc = centroid(Q)
    P -= Pc
    Q -= Qc

    # calculate the covrariance matrix A
    A = np.dot(np.transpose(P), Q)

    # calculate the SVD components
    V, S, W = np.linalg.svd(A)

    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # calculate rotation matrix U
    U = np.dot(V, W)

    return A, V, S, W, d, U

def kabsch_rotate(P, Q):
    """
    Rotate matrix P unto matrix Q using kabsch algorithm
    :param P:
    :param Q:
    :return:
    """
    U = kabsch(P, Q)[-1]
    P_rot = np.dot(P, U)
    return P_rot

def kabsch_rmsd(P, Q):
    """
    Rotate the matrix P unto Q using function kabsch_rotate and return the rmsd calculated
    :param P:
    :param Q:
    :return:
    """
    P = kabsch_rotate(P, Q)
    rmsd_ = rmsdcalc(P, Q)
    return rmsd_

def calculate_rmsd(file1, file2):
    file1_coord = get_coordinates(file1)
    file2_coord = get_coordinates(file2)
    rmsd_ = kabsch_rmsd(file1_coord, file2_coord)
    return rmsd_

if __name__ == '__main__':
    file1 = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\ab40_c0_ySERF_01_2646.pdb"
    file2 = r"C:\Users\sugyan\Documents\MembraneMDfiles\serf_ab_models\ab40_c0_ySERF_01_5677.pdb"
    rmsd = calculate_rmsd(file1, file2)
    print(rmsd)