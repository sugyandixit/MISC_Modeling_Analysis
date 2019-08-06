# conv_cg_aa_to_pdb.py

import os


def files(dirpath):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith('cgaa.pdb')]
    return files


def writepdb(fname, list, mode):
    """
    Writes pdb file using input list that contains all the field inputs for pdb file

    Parameters
    ------
    fname: string
        .pdb file to read
    list: array [n_num_atoms, n_fields]
    n_fields = atm, atm_sernum, atm_altloc, atm_resname, atm_chainid, atm_resseqnum, atm_cod, atm_xcoor, atm_ycoor,
    atm_zcoor, atm_occ, atm_tempf, atm_segid, atm_elem, atm_charge

    Returns
    -------
    self: object

    """

    file = open(fname, mode='a' or 'w')

    for i in range(len(list[0])):
        atm = list[0][i].ljust(6)
        # print (atm)
        atm_sernum = list[1][i].rjust(5)
        atm_space0 = ' '
        atm_name = list[2][i].ljust(4)
        atm_altloc = list[3][i]
        atm_resname = list[4][i].rjust(3)
        atm_space1 = ' '
        atm_chainid = list[5][i]
        atm_resseqnum = list[6][i].rjust(4)
        atm_cod = list[7][i]
        atm_space2 = '   '
        atm_xcoor = str('%8.3f' % (float(list[8][i]))).rjust(8)
        atm_ycoor = str('%8.3f' % (float(list[9][i]))).rjust(8)
        atm_zcoor = str('%8.3f' % (float(list[10][i]))).rjust(8)
        atm_occ = str('%6.2f' % (float(list[11][i]))).rjust(6)
        atm_tempf = str('%6.2f' % (float(list[12][i]))).rjust(6)
        atm_space3 = '      '
        atm_segid = list[13][i].rjust(4)
        atm_elem = list[14][i].rjust(2)
        atm_charge = list[15][i].ljust(2)

        # print (atm_xcoor)

        # file.write("%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n" % atm, atm_sernum, atm_space0, atm_name, atm_altloc,
        # atm_resname, atm_space1, atm_chainid, atm_resseqnum, atm_cod, atm_space2, atm_xcoor, atm_ycoor, atm_zcoor,
        # atm_occ, atm_tempf, atm_space3, atm_segid, atm_elem, atm_charge)

        file.write(
            atm + atm_sernum + atm_space0 + atm_name + atm_altloc + atm_resname + atm_space1 + atm_chainid +
            atm_resseqnum + atm_cod + atm_space2 + atm_xcoor + atm_ycoor + atm_zcoor + atm_occ + atm_tempf +
            atm_space3 + atm_segid + atm_elem + atm_charge + '\n')

    ter = 'TER   ' + str(int(list[1][-1]) + 1).rjust(5) + '      ' + str(list[4][-1]).rjust(3) + ' ' + str(
        list[5][-1]) + list[6][-1].rjust(4) + ' \n'
    file.write(ter)

    file.close()

def write_xyz_file(fname, atomsymb, xcoords, ycoords, zcoords):

    with open(fname, 'w') as xyzfile:
        header = 'xyz file\n'
        xyzfile.write(header)
        for ind, atom in enumerate(atomsymb):
            line = '{}   {}   {}   {}\n'.format(atom, xcoords[ind], ycoords[ind], zcoords[ind])
            xyzfile.write(line)



if __name__=='__main__':
    dirpath = r"C:\Users\sugyan\Desktop"
    fs = files(dirpath)
    for file in fs:
        fname_0 = str(file).split('.')[0]
        fname = str(file).split('.')[0] + '_pdbout.pdb'
        print(fname)
        cgaa = open(os.path.join(dirpath, file), 'r+')
        cgaa = cgaa.read().splitlines()
        atm = []
        atm_sernum = []
        atm_name = []
        atm_altloc = []
        atm_resname = []
        atm_chainid = []
        atm_resseqnum = []
        atm_cod = []
        atm_xcoor = []
        atm_ycoor = []
        atm_zcoor = []
        atm_occ = []
        atm_tempf = []
        atm_segid = []
        atm_elem = []
        atm_charge = []
        for line in cgaa[2:-1]:
            atm.append('ATOM')
            atm_sernum1 = line[15:20].strip(' ')
            atm_sernum.append(atm_sernum1)
            atm_name1 = line[11:15].strip(' ')
            atm_name.append(atm_name1)
            atm_altloc.append(' ')
            atm_resname1 = line[5:9].strip(' ')
            atm_resname.append(atm_resname1)
            atm_chainid.append(' ')
            atm_resseqnum1 = line[1:5].strip(' ')
            atm_resseqnum.append(atm_resseqnum1)
            atm_cod.append(' ')
            xcoor = line[20:28].strip(' ')
            ycoor = line[28:36].strip(' ')
            zcoor = line[36:44].strip(' ')
            atm_xcoor.append(xcoor)
            atm_ycoor.append(ycoor)
            atm_zcoor.append(zcoor)
            atm_occ.append(1.0)
            atm_tempf.append(0.0)
            atm_segid.append(' ')
            atm_elem.append(atm_name1[0])
            atm_charge.append(' ')

        list = [atm, atm_sernum, atm_name, atm_altloc, atm_resname, atm_chainid, atm_resseqnum, atm_cod, atm_xcoor,
                atm_ycoor, atm_zcoor, atm_occ, atm_tempf, atm_segid, atm_elem, atm_charge]

        # writepdb(os.path.join(dirpath, fname), list, 'w')

        xyzfname = fname_0+'_xyzfile.xyz'
        write_xyz_file(os.path.join(dirpath, xyzfname), atm_elem, atm_xcoor, atm_ycoor, atm_zcoor)