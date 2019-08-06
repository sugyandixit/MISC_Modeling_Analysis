# convertcharmpdbtoregpdb.py

import os


def inputfiles(dirpath):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('.pdb')]
    return files


def conv_pdb_for_imos(dirpath):
    files = inputfiles(dirpath)

    outputdirname = '/convcharmmpdbtoregpdb_output'
    outputpath = os.path.join(dirpath + outputdirname)
    if not os.path.isdir(outputpath): os.makedirs(outputpath)

    for file in files:
        outfilename1 = str(file).split('.')[0] + '_out.pdb'
        outfilename = os.path.join(outputpath, outfilename1)
        outputfile = open(outfilename, 'w')
        charmmpdb = open(os.path.join(dirpath, file), 'r+')
        charmmpdb = charmmpdb.read().splitlines()
        data = []
        atomid1 = []
        # outputfile = open("out.pdb", 'w')
        x = []
        for line in charmmpdb:
            # print(len(line))
            if line.startswith('ATOM'):
                letter = line[13]
                if letter == 'D':
                    letter = 'H'
                if letter == 'E':
                    letter = 'H'
                line2 = line + " " + letter + '\n'
                outputfile.write(line2)
            if line.startswith('TER'):
                line3 = line + '\n'
                outputfile.write(line3)
            if line.startswith('END'):
                outputfile.write(line)
        outputfile.close()


if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\IMOS\models\serf_nomod_structs_random_selection"
    conv_pdb_for_imos(dirpath)