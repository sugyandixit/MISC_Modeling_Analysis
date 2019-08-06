# imos_setup_v1.py
import os
import subprocess
import time


def make_new_dir(dirpath):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)
    return dirpath

def readmodels(dirpath):
    models = os.listdir(dirpath)
    models = [x for x in models if x.endswith('.pdb')]
    return models


def setupIMOSfile(IMOS_dirpath, input_pdb_fpath, output_fpath):
    filename = 'IMoS.cla'
    imosfileread = open(os.path.join(IMOS_dirpath, filename), 'r')
    imfile = imosfileread.read().splitlines()
    imfile[1] = str(input_pdb_fpath) + ' ' + str(output_fpath) + '  He'
    mass, charge = calcmassandcharge(input_pdb_fpath)
    imfile[7] = 'Charge ' + str(charge)
    imfile[12] = 'Mweight ' + str(mass)
    imosfile_write = open(os.path.join(IMOS_dirpath, filename), 'w')
    imosfile_write.write('\n'.join(imfile))


def calcmassandcharge(input_pdb_fpath):
    mass_list = [['C', 'N', 'O', 'H', 'F', 'P', 'S'],
                 [12.0107, 14.0067, 15.9994, 1.00794, 18.998403, 30.973762, 32.065]]

    charmmpdb = open(input_pdb_fpath, 'r+')
    charmmpdb = charmmpdb.read().splitlines()
    mass = 0
    charge = 0
    for line in charmmpdb:
        if line.startswith('ATOM'):
            atom = line[77]
            for i in range(len(mass_list[0])):
                if atom == mass_list[0][i]:
                    atommass = float(mass_list[1][i])
                    mass = atommass + mass
            partial_charge = float(line[61:66])
            charge += partial_charge
    charge = int(charge)
    return round(mass, 4), charge


def run_imos(IMOS_dirpath, input_pdb_dirpath, output_dir):

    start_time_init = time.perf_counter()

    input_pdb_file_list = readmodels(input_pdb_dirpath)
    output_dirname = os.path.join(os.getcwd(), output_dir)
    output_dirpath = make_new_dir(output_dirname)


    for index, pdb_name in enumerate(input_pdb_file_list):
        st_time_2 = time.perf_counter()

        input_pdb_fpath = os.path.join(input_pdb_dirpath, pdb_name)

        output_fname = pdb_name + '_out.txt'
        output_fpath = os.path.join(output_dir, output_fname)

        imos_file = setupIMOSfile(IMOS_dirpath, input_pdb_fpath, output_fpath)
        subprocess.call(os.path.join(IMOS_dirpath, "IMoS106W64d.exe"))

        time_now = time.perf_counter()

        time_took = time_now - st_time_2

        print("---------- %s seconds gone----------" % time_took)

        st_time_2 = time.perf_counter()

    total_time = time.perf_counter() - start_time_init

    print("---------- %s seconds Total ----------" % (time.clock() - start_time_init))




if __name__=='__main__':
    imos_dirpath = os.getcwd()
    input_pdb_dirpath = 'models\serf_mod_structs_random_selection\convcharmmpdbtoregpdb_output'
    output_dirname = 'savefolder\serf_mod_structs_random_selection'
    run_imos(imos_dirpath, input_pdb_dirpath, output_dirname)



