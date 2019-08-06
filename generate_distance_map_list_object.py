import os
import pickle
import distance_map
from distance_map import generate_distance_map, plot_contour_distmat
# from dataclasses import dataclass

# @dataclass
class DistanceMap(object):
    """
    DistanceMap Object to store info on dirpath, pdbname, chain_ids, and distance matrix
    """

    def __init__(self):
        """
        Initialize
        """

        self.dirpath = None
        self.pdbname = None
        self.chain_id = None
        self.distance_mat = None


def get_pdb_files(dirpath):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('.pdb')]
    return files


def save_object_to_pickle(obj, file_path):
    with open(file_path, 'wb') as out_obj:
        pickle.dump(obj, out_obj)
        out_obj.close()


if __name__ == '__main__':
    dirpath = r'C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\step2_remd_coor_final'
    pdbfiles = get_pdb_files(dirpath)
    chain_ids = None
    dirpath_list = []
    pdb_list = []
    chainid_list = []
    distmap_list = []
    for pdb in pdbfiles:
        dist_map = generate_distance_map(pdb, dirpath, chain_ids)
        plot_contour_distmat(dist_map, pdb, dirpath)
        dirpath_list.append(dirpath)
        pdb_list.append(pdb)
        chainid_list.append(chain_ids)
        distmap_list.append(dist_map)
    dist_map_obj = DistanceMap()
    dist_map_obj.dirpath = dirpath_list
    dist_map_obj.pdbname = pdb_list
    dist_map_obj.chain_id = chainid_list
    dist_map_obj.distance_mat = distmap_list

    out_object_path = os.path.join(dirpath, 'distance_matrix.obj')
    save_object_to_pickle(dist_map_obj, out_object_path)
    print('heho')