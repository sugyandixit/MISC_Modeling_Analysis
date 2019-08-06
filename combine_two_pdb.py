from pdb_utility import *


directory = r"C:\Users\sugyan\Documents\MembraneMDfiles\PNP_dimer_membprot"
pdb1 = 'PMPwSS.pdb'
pdb2 = 'PMPwSS.pdb'

pdb1_obj = ParsePDB(pdb1, directory).read()
pdb2_obj = ParsePDB(pdb2, directory).read()

#translate and/or rotate coordinates from pdb1
xcoor1 = np.array(pdb1_obj.atm_xcoor, dtype='float')
ycoor1 = np.array(pdb1_obj.atm_ycoor, dtype='float')
zcoor1 = np.array(pdb1_obj.atm_zcoor, dtype='float')

com_1 = centerofmass(xcoor1, ycoor1, zcoor1)

coors_1 = mergecoords(xcoor1, ycoor1, zcoor1)
coors_1_trans = translate(coors_1, delta_x=30, delta_z=-10, delta_y=20)
# rot_matrix = rotation_matrix(10, direction='x')
# coors_1_rotate = rotate_structure(coors_1_trans, rot_matrix)
coors_1_final = coors_1_trans

com_1_final = centerofmass(coors_1_final[:, 0], coors_1_final[:, 1], coors_1_final[:, 2])

#translate and/or rotate coordinates from pdb2
xcoor2 = np.array(pdb2_obj.atm_xcoor, dtype='float')
ycoor2 = np.array(pdb2_obj.atm_ycoor, dtype='float')
zcoor2 = np.array(pdb2_obj.atm_zcoor, dtype='float')

com_2 = centerofmass(xcoor2, ycoor2, zcoor2)


com_distance = distance_between_two_points(com_1, com_2)
print(com_distance)


coors_2 = mergecoords(xcoor2, ycoor2, zcoor2)
coors_2_trans = translate(coors_2, delta_x=-1)
coors_2_final = coors_2_trans

com_2_final = centerofmass(coors_2_final[:, 0], coors_2_final[:, 1], coors_2_final[:, 2])

com_distance_final = distance_between_two_points(com_1_final, com_2_final)
print(com_distance_final)

coor_list = [coors_1_final, coors_2_final]
pdb_list = [pdb1_obj, pdb2_obj]
chainid_list = ['A', 'S']

list_of_fields = []
i=1
for index, (pdb, chainid, coors) in enumerate(zip(pdb_list, chainid_list, coor_list)):
    xcoor = coors[:, 0]
    ycoor = coors[:, 1]
    zcoor = coors[:, 2]
    for ind, atm in enumerate(pdb.atm):
        list_of_fields.append([atm, str(i), pdb.atm_name[ind], pdb.atm_altloc[ind], pdb.atm_resname[ind],
                               chainid, pdb.atm_resseqnum[ind], pdb.atm_cod[ind], xcoor[ind],
                               ycoor[ind], zcoor[ind], pdb.atm_occ[ind], pdb.atm_tempf[ind],
                               pdb.atm_segid[ind], pdb.atm_elem[ind], pdb.atm_charge[ind]])
        i += 1


fname_dock = pdb1.split('.')[0] + '_' + pdb2.split('.')[0] + '_dock_init.pdb'

write_pdb(list_of_fields, fname_dock, directory, elem_charge=True)

print('heho')