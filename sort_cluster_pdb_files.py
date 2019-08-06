# put cluster pdbs in cluster directories

import os
import subprocess

def make_directory(dirpath):

	if not os.path.exists(dirpath):
		os.makedirs(dirpath)
	
	return dirpath


def sort_cluster_pdb_files(dirpath, fname, out_cluster_dirname):
	
	cluster_num = []
	file_path = []
	
	with open(os.path.join(dirpath, fname), 'r') as clust_pdb_file:
		
		fileread = clust_pdb_file.read().splitlines()
		for line in fileread[1:]:
			chars = line.split(',')
			print('cluster ', chars[0])
			print('pdb_file', chars[1])
			clust_dir_name = out_cluster_dirname + '_' + str(chars[0])
			clust_dir_path = make_directory(os.path.join(dirpath, clust_dir_name))
			copy_command_args = ['cp', chars[1], str(clust_dir_path)+'/']
			subprocess.call(copy_command_args)
			
if __name__ == '__main__':

	dirpath = os.getcwd()
	fname_auto_cluster = 'cluster_auto_pdb_files.csv'
	fname_manual_cluster = 'cluster_manual_pdb_files.csv'
	#sort_cluster_pdb_files(dirpath, fname_auto_cluster, 'cluster_auto')
	sort_cluster_pdb_files(dirpath, fname_manual_cluster, 'cluster_manual')