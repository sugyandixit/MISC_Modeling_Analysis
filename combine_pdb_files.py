# combine pdbs in a directory to a single pdb file

import os

def get_pdb_files(dirpath, startid):
	
	files = os.listdir(dirpath)
	files = [x for x in files if x.endswith('.pdb')]
	files = [x for x in files if x.startswith(startid)]
	return files
	

def combine_pdb_files(list_of_pdb_files, dirpath):
	
	output_string = ''
	header = 'REMARK  COMBINED  PDB  MODELS\n'
	output_string += header
	
	for ind, file in enumerate(list_of_pdb_files):
		
		model_header = 'MODEL '+ str(ind+1) + '\n'
		output_string += model_header
		
		with open(os.path.join(dirpath, file), 'r') as pdbfile:
			pdb_read = pdbfile.read().splitlines()
			for line in pdb_read:
				if line.startswith('ATOM'):
					output_string += line + '\n'
		end_mod_header = 'TER \n'
		end_model = 'ENDMDL \n'
		output_string += end_mod_header
		output_string += end_model
	
	
	combo_pdb_outname = 'combo_pdb.pdb'
	
	with open(os.path.join(dirpath, combo_pdb_outname), 'w') as outfile:
		outfile.write(output_string)

if __name__ == '__main__':
	
	# dirpath = os.getcwd()
	dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\charmm_gui_iapp_nanodisc\iappdimer_integral\19\step7_structs\output\combine_pdb_files"
	startid = 'step7'
	list_pdb_files = get_pdb_files(dirpath, startid)
	combine_pdb_files(list_pdb_files, dirpath)
