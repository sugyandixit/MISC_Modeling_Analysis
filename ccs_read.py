import os

def get_impccs_files(dirpath):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('.impccs')]
    return files

def read_impact_ccs_file(file, dirpath):
    ccsoutfile = open(os.path.join(dirpath, file), 'r+')
    ccsoutfile = ccsoutfile.read().splitlines()

    fname = []
    CCS_PA = []
    CCS_TM = []

    for line in ccsoutfile:
        if line.startswith('CCS'):
            stt, num, filename, ccspa, obj1, semrel1, ccstm, niter = line.split()
            ccspa = float(ccspa)
            ccstm = float(ccstm)
            CCS_PA.append(ccspa)
            CCS_TM.append(ccstm)
            fname.append(filename)

    return fname, CCS_PA, CCS_TM



def conv_imp_to_csv(impccsfname, imp_fname, ccs_pa, ccs_tm, dirpath, frame=None):
    impccsfname_csv = str(impccsfname).split('.impccs')[0] + '.csv'
    with open(os.path.join(dirpath, impccsfname_csv), 'w') as impccs:
        if frame==None:
            header = '{},{},{},{}\n'.format('number','filename','CCS_PA','CCS_TM')
            impccs.write(header)
            for ind in range(len(imp_fname)):
                line = '{},{},{},{}\n'.format(ind, imp_fname[ind], ccs_pa[ind], ccs_tm[ind])
                impccs.write(line)
        else:
            header = '{},{},{},{},{}\n'.format('number','frame','filename','CCS_PA','CCS_TM')
            impccs.write(header)
            for ind in range(len(imp_fname)):
                line = '{},{},{},{},{}\n'.format(ind, frame[ind], imp_fname[ind], ccs_pa[ind], ccs_tm[ind])
                impccs.write(line)
    impccs.close()



def read_impccs_to_csv(impfname, dirpath, frame=None):
    fname, ccs_pa, ccs_tm = read_impact_ccs_file(impfname, dirpath)
    if frame:
        frame_ = []
        for item in fname:
            chars = item.split('.pdb')[0]
            framenum = chars.split('_')[-1]
            frame_.append(framenum)
        conv_imp_to_csv(impfname, fname, ccs_pa, ccs_tm, dirpath, frame=frame_)
    else:
        conv_imp_to_csv(impfname, fname, ccs_pa, ccs_tm, dirpath, frame=frame)

if __name__=='__main__':
    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\PNP_dimer_membprot"
    impccs_files = get_impccs_files(dirpath)
    for file in impccs_files:
        read_impccs_to_csv(file, dirpath, frame=None)