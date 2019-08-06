import os
import numpy as np
from read_dssp_file import DSSPObject
from generate_dssp_list_object import DSSPListObj
from generate_dssp_list_object import load_pickle_object, save_object_to_pickle


def generate_secondary_struct_dict():
    """
    generate secondary structure dictionary with values corresponding to output from dssp
    :return: nested dictionary of secondary structure [helix, beta_strand, loop]
    """

    helix_dict = dict()
    helix_dict['310_helix'] = 'G'
    helix_dict['alpha_helix'] = 'H'
    helix_dict['pi_helix'] = 'I'

    beta_dict = dict()
    beta_dict['beta_buldge'] = 'E'
    beta_dict['beta_bridge'] = 'B'

    loop_dict = dict()
    loop_dict['turns'] = 'T'
    loop_dict['curvature'] = 'S'
    loop_dict['other_loop'] = ' '.isspace()

    sec_struct_dict = dict()
    sec_struct_dict['helix'] = helix_dict
    sec_struct_dict['beta_strand'] = beta_dict
    sec_struct_dict['loop'] = loop_dict

    return sec_struct_dict


def dssp_obj_analysis_dictionary(frac_helix, frac_beta_strand, frac_loop, frac_310_helix, frac_alpha_helix, frac_pi_helix,
                                 frac_beta_buldge, frac_beta_bridge, frac_turns, frac_curve, frac_other_loop, bridge_res,
                                      bridge_partner_1, bridge_partner_2, bridge_label):
    """
    dictionary to store fraction sec struct information along with bridge partner info
    :param frac_helix: fraction helix including all types
    :param frac_beta_strand: fraction beta strand including all types
    :param frac_loop: fraction loop including all types
    :param frac_310_helix: fraction 310 helix
    :param frac_alpha_helx: fraction alpha helix
    :param frac_pi_helix: fraction pi helix
    :param frac_beta_buldge: fraction beta buldge
    :param frac_beta_bridge: fraction beta bridge
    :param frac_turns: fraction turns
    :param frac_curve: fraction curves
    :param frac_other_loop: fraction other loops
    :param brdige_res: array of bridge res
    :param bridge_partner_1: array of bridge partners 1
    :param bridge_partner_2: array of bridge partners 2
    :param bridge_label: array of bridge labels
    :return: dictionary
    """
    dssp_analys_dict = dict()
    dssp_analys_dict['frac_helix'] = frac_helix
    dssp_analys_dict['frac_beta_strand'] = frac_beta_strand
    dssp_analys_dict['frac_loop'] = frac_loop
    dssp_analys_dict['frac_310_helix'] = frac_310_helix
    dssp_analys_dict['frac_alpha_helix'] = frac_alpha_helix
    dssp_analys_dict['frac_pi_helix'] = frac_pi_helix
    dssp_analys_dict['frac_beta_buldge'] = frac_beta_buldge
    dssp_analys_dict['frac_beta_bridge'] = frac_beta_bridge
    dssp_analys_dict['frac_turns'] = frac_turns
    dssp_analys_dict['frac_curve'] = frac_curve
    dssp_analys_dict['frac_other_loop'] = frac_other_loop
    dssp_analys_dict['bridge_res'] = bridge_res
    dssp_analys_dict['bridge_partner_1'] = bridge_partner_1
    dssp_analys_dict['bridge_partner_2'] = bridge_partner_2
    dssp_analys_dict['bridge_label'] = bridge_label

    return dssp_analys_dict


def dssp_obj_analysis_struct(dssp_obj):
    """
    determine fraction sec structues, num_hbonds, accessible surface area, and bridge partners from dssp object
    :param dssp_obj: dssp object
    :return: dssp_analysis_dictionary
    """

    len_of_res = len(dssp_obj.struct_code)

    num_310_helix = 0
    num_alpha_helix = 0
    num_pi_helix = 0
    num_beta_buldge = 0
    num_beta_bridge = 0
    num_turns = 0
    num_curve = 0
    num_other_loop = 0

    bridge_res = []
    bridge_partner_1 = []
    bridge_partner_2 = []
    bridge_label = []

    sec_struct_dict = generate_secondary_struct_dict()
    helix_dict = sec_struct_dict['helix']
    beta_strand_dict = sec_struct_dict['beta_strand']
    loop_dict = sec_struct_dict['loop']

    for ind, (sec_struc_code, bridge_lab) in enumerate(zip(dssp_obj.struct_code, dssp_obj.bridge_label)):

        if sec_struc_code == helix_dict['310_helix']:
            num_310_helix += 1

        if sec_struc_code == helix_dict['alpha_helix']:
            num_alpha_helix += 1

        if sec_struc_code == helix_dict['pi_helix']:
            num_pi_helix += 1

        if sec_struc_code == beta_strand_dict['beta_buldge']:
            num_beta_buldge += 1

        if sec_struc_code == beta_strand_dict['beta_bridge']:
            num_beta_bridge += 1

        if sec_struc_code == loop_dict['turns']:
            num_turns += 1

        if sec_struc_code == loop_dict['curvature']:
            num_curve += 1

        if sec_struc_code.isspace() == loop_dict['other_loop']:
            num_other_loop += 1

        if not bridge_lab.isspace():
            bridge_res.append(str(dssp_obj.res_num[ind])+ '_' + dssp_obj.aa_code[ind])
            bridge_label.append(bridge_lab)
            if dssp_obj.bridge_partner_1[ind] > 0:
                bridge_partner_1.append(str(dssp_obj.bridge_partner_1[ind]) + '_' + dssp_obj.aa_code[dssp_obj.bridge_partner_1[ind]-1])
            else:
                bridge_partner_1.append(dssp_obj.bridge_partner_1[ind])
            if dssp_obj.bridge_partner_2[ind] > 0:
                bridge_partner_2.append(str(dssp_obj.bridge_partner_2[ind]) + '_' + dssp_obj.aa_code[dssp_obj.bridge_partner_2[ind]-1])
            else:
                bridge_partner_2.append(dssp_obj.bridge_partner_2[ind])

    total_num_helix = num_310_helix + num_alpha_helix + num_pi_helix
    total_num_beta_strand = num_beta_buldge + num_beta_bridge
    total_num_loop = num_turns + num_curve + num_other_loop


    dssp_analys_dict = dssp_obj_analysis_dictionary(frac_helix=total_num_helix/len_of_res,
                                                    frac_beta_strand =total_num_beta_strand/len_of_res,
                                                    frac_loop=total_num_loop/len_of_res,
                                                    frac_310_helix=num_310_helix/len_of_res,
                                                    frac_alpha_helix=num_alpha_helix/len_of_res,
                                                    frac_pi_helix=num_pi_helix/len_of_res,
                                                    frac_beta_buldge=num_beta_buldge/len_of_res,
                                                    frac_beta_bridge=num_beta_bridge/len_of_res,
                                                    frac_turns=num_turns/len_of_res,
                                                    frac_curve=num_curve/len_of_res,
                                                    frac_other_loop=num_other_loop/len_of_res,
                                                    bridge_res=bridge_res,
                                                    bridge_partner_1=bridge_partner_1,
                                                    bridge_partner_2=bridge_partner_2,
                                                    bridge_label=bridge_label)

    return dssp_analys_dict


def write_dssp_analysis_struct(dssp_list_obj, dirpath, outid):
    """
    write dssp analysis to a file. For one dssp object. Another writing module for list of dssp objects
    :param dssp_obj: single dssp object
    :param dssp_analysis_dict: dssp analysis dict
    :param dirpath: directory to save file
    :return:void. Writes and saves file
    """

    header_string = 'dssp_fname,accessible_surface_area,total_num_hbonds,hbonds_parallel_bridge,hbonds_antiparallel_bridge,hbonds_0plus,hbonds_i1_minus,hbonds_i1_plus,hbonds_i2_minus,hbonds_i2_plus,hbonds_i3_minus,hbonds_i3_plus,hbonds_i4_minus,hbonds_i4_plus,hbonds_i5_minus,hbonds_i5_plus,frac_helix,frac_beta_strand,frac_loop,frac_310_helix,frac_alpha_helix,frac_pi_helix,frac_beta_buldge,frac_beta_bridge,frac_turns,frac_curve,frac_other_loop\n'

    data_string = ''

    for index, dssp_obj in enumerate(dssp_list_obj.dssp_object):

        dssp_analysis_dict = dssp_obj_analysis_struct(dssp_obj)
        dssp_fname = os.path.split(dssp_obj.dssp_fpath)[1]

        line = '{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(
            dssp_fname,
            dssp_obj.accessible_surface,
            dssp_obj.total_num_hbonds,
            dssp_obj.total_num_hbonds_parallel_bridge,
            dssp_obj.total_num_hbonds_antiparallel_bridge,
            dssp_obj.total_num_hbonds_i0_plus,
            dssp_obj.total_num_hbonds_i1_minus,
            dssp_obj.total_num_hbonds_i1_plus,
            dssp_obj.total_num_hbonds_i2_minus,
            dssp_obj.total_num_hbonds_i2_plus,
            dssp_obj.total_num_hbonds_i3_minus,
            dssp_obj.total_num_hbonds_i3_plus,
            dssp_obj.total_num_hbonds_i4_minus,
            dssp_obj.total_num_hbonds_i4_plus,
            dssp_obj.total_num_hbonds_i5_minus,
            dssp_obj.total_num_hbonds_i5_plus,
            dssp_analysis_dict['frac_helix'],
            dssp_analysis_dict['frac_beta_strand'],
            dssp_analysis_dict['frac_loop'],
            dssp_analysis_dict['frac_310_helix'],
            dssp_analysis_dict['frac_alpha_helix'],
            dssp_analysis_dict['frac_pi_helix'],
            dssp_analysis_dict['frac_beta_buldge'],
            dssp_analysis_dict['frac_beta_bridge'],
            dssp_analysis_dict['frac_turns'],
            dssp_analysis_dict['frac_curve'],
            dssp_analysis_dict['frac_other_loop'])#,
            # dssp_analysis_dict['bridge_res'],
            # dssp_analysis_dict['bridge_partner_1'],
            # dssp_analysis_dict['bridge_partner_2'],
            # dssp_analysis_dict['bridge_label'])

        data_string += line

    output_string = header_string + data_string


    out_name = os.path.join(dirpath, outid + '_dssp_analysis_struct.csv')

    with open(out_name, 'w') as outfile:
        outfile.write(output_string)
        outfile.close()


def dssp_list_obj_analysis_dictionary(res_array, frac_helix, frac_beta_strand, frac_loop, frac_310_helix, frac_alpha_helix, frac_pi_helix,
                                 frac_beta_buldge, frac_beta_bridge, frac_turns, frac_curve, frac_other_loop):
    """
    dictionary to store fraction sec struct information along with bridge partner info
    :param res_array: residue array
    :param frac_helix: fraction helix including all types
    :param frac_beta_strand: fraction beta strand including all types
    :param frac_loop: fraction loop including all types
    :param frac_310_helix: fraction 310 helix
    :param frac_alpha_helx: fraction alpha helix
    :param frac_pi_helix: fraction pi helix
    :param frac_beta_buldge: fraction beta buldge
    :param frac_beta_bridge: fraction beta bridge
    :param frac_turns: fraction turns
    :param frac_curve: fraction curves
    :param frac_other_loop: fraction other loops
    :return: dictionary
    """
    dssp_analys_dict = dict()
    dssp_analys_dict['res_array'] = res_array
    dssp_analys_dict['frac_helix'] = frac_helix
    dssp_analys_dict['frac_beta_strand'] = frac_beta_strand
    dssp_analys_dict['frac_loop'] = frac_loop
    dssp_analys_dict['frac_310_helix'] = frac_310_helix
    dssp_analys_dict['frac_alpha_helix'] = frac_alpha_helix
    dssp_analys_dict['frac_pi_helix'] = frac_pi_helix
    dssp_analys_dict['frac_beta_buldge'] = frac_beta_buldge
    dssp_analys_dict['frac_beta_bridge'] = frac_beta_bridge
    dssp_analys_dict['frac_turns'] = frac_turns
    dssp_analys_dict['frac_curve'] = frac_curve
    dssp_analys_dict['frac_other_loop'] = frac_other_loop

    return dssp_analys_dict



def write_dssp_analysis_res(dssp_list_obj, dirpath, outid):
    """
    write the sec structure frac for each residue
    :param dssp_list_obj: dssp_list_obj
    :return: void. Writes file
    """

    dssp_analysis_dict = dssp_list_obj_analysis_res(dssp_list_obj)

    header = 'res_code,res_num,frac_helix,frac_beta_strand,frac_loop,frac_310_helix,frac_alpha_helix,frac_pi_helix,frac_beta_buldge,frac_beta_bridge,frac_turns,frac_curve,frac_other_loops\n'

    data_string = ''

    for num in range(len(dssp_analysis_dict['res_array'])):
        line = '{},{},{},{},{},{},{},{},{},{},{},{},{}\n'.format(dssp_analysis_dict['res_array'][num],
                                                                 num + 1,
                                                                 dssp_analysis_dict['frac_helix'][num],
                                                                 dssp_analysis_dict['frac_beta_strand'][num],
                                                                 dssp_analysis_dict['frac_loop'][num],
                                                                 dssp_analysis_dict['frac_310_helix'][num],
                                                                 dssp_analysis_dict['frac_alpha_helix'][num],
                                                                 dssp_analysis_dict['frac_pi_helix'][num],
                                                                 dssp_analysis_dict['frac_beta_buldge'][num],
                                                                 dssp_analysis_dict['frac_beta_bridge'][num],
                                                                 dssp_analysis_dict['frac_turns'][num],
                                                                 dssp_analysis_dict['frac_curve'][num],
                                                                 dssp_analysis_dict['frac_other_loop'][num])
        data_string += line


    output_string = header + data_string

    out_name = os.path.join(dirpath, outid + '_dssp_analysis_res.csv')

    with open(out_name, 'w') as outfile:
        outfile.write(output_string)
        outfile.close()





def dssp_list_obj_analysis_res(dssp_list_obj):
    """
    take all dssp object in list to determine sec structure for each residue
    :param dssp_list_obj: dssp list object
    :return: dictionary
    """

    sec_struct_array = []

    res_array = dssp_list_obj.dssp_object[0].aa_code

    for index, dssp_obj in enumerate(dssp_list_obj.dssp_object):

        sec_struct = dssp_obj.struct_code
        sec_struct_array.append(sec_struct)




    sec_struct_res_array = np.array(sec_struct_array).T

    sec_struct_dict = generate_secondary_struct_dict()
    helix_dict = sec_struct_dict['helix']
    beta_strand_dict = sec_struct_dict['beta_strand']
    loop_dict = sec_struct_dict['loop']

    frac_helix_res = []
    frac_beta_strand_res = []
    frac_loop_res = []
    frac_310_helix_res = []
    frac_alpha_helix_res = []
    frac_pi_helix_res = []
    frac_beta_buldge_res = []
    frac_beta_bridge_res = []
    frac_turns_res = []
    frac_curve_res = []
    frac_other_loop_res = []


    for index, res_struct_arr in enumerate(sec_struct_res_array):

        num_310_helix = 0
        num_alpha_helix = 0
        num_pi_helix = 0
        num_beta_buldge = 0
        num_beta_bridge = 0
        num_turns = 0
        num_curve = 0
        num_other_loop = 0

        for ind, sec_struc_code in enumerate(res_struct_arr):

            if sec_struc_code == helix_dict['310_helix']:
                num_310_helix += 1

            if sec_struc_code == helix_dict['alpha_helix']:
                num_alpha_helix += 1

            if sec_struc_code == helix_dict['pi_helix']:
                num_pi_helix += 1

            if sec_struc_code == beta_strand_dict['beta_buldge']:
                num_beta_buldge += 1

            if sec_struc_code == beta_strand_dict['beta_bridge']:
                num_beta_bridge += 1

            if sec_struc_code == loop_dict['turns']:
                num_turns += 1

            if sec_struc_code == loop_dict['curvature']:
                num_curve += 1

            if sec_struc_code.isspace() == loop_dict['other_loop']:
                num_other_loop += 1

        total_num_helix = num_310_helix + num_alpha_helix + num_pi_helix
        total_num_beta_strand = num_beta_buldge + num_beta_bridge
        total_num_loop = num_turns + num_curve + num_other_loop

        len_structs = len(res_struct_arr)

        frac_helix_res.append(total_num_helix/len_structs)
        frac_beta_strand_res.append(total_num_beta_strand/len_structs)
        frac_loop_res.append(total_num_loop/len_structs)
        frac_310_helix_res.append(num_310_helix/len_structs)
        frac_alpha_helix_res.append(num_alpha_helix/len_structs)
        frac_pi_helix_res.append(num_pi_helix/len_structs)
        frac_beta_buldge_res.append(num_beta_buldge/len_structs)
        frac_beta_bridge_res.append(num_beta_bridge/len_structs)
        frac_turns_res.append(num_turns/len_structs)
        frac_curve_res.append(num_curve/len_structs)
        frac_other_loop_res.append(num_other_loop/len_structs)

    analysis_dict = dssp_list_obj_analysis_dictionary(res_array=res_array,
                                                      frac_helix=frac_helix_res,
                                                      frac_beta_strand=frac_beta_strand_res,
                                                      frac_loop=frac_loop_res,
                                                      frac_310_helix=frac_310_helix_res,
                                                      frac_alpha_helix=frac_alpha_helix_res,
                                                      frac_pi_helix=frac_pi_helix_res,
                                                      frac_beta_buldge=frac_beta_buldge_res,
                                                      frac_beta_bridge=frac_beta_bridge_res,
                                                      frac_turns=frac_turns_res,
                                                      frac_curve=frac_curve_res,
                                                      frac_other_loop=frac_other_loop_res)

    print('heho')

    return analysis_dict


def get_dssp_list_object_files(dirpath, endid):
    files = os.listdir(dirpath)
    files = [x for x in files if x.endswith(endid)]
    return files

def dssp_list_object_analysis_all(dirpath, endid='_dssp_list_object.obj'):

    dssp_list_obj_files = get_dssp_list_object_files(dirpath, endid)

    for ind, file in enumerate(dssp_list_obj_files):

        outid = str(file).split(endid)[0]

        dssp_list_obj = load_pickle_object(os.path.join(dirpath, file))

        write_dssp_analysis_res(dssp_list_obj, dirpath, outid)
        write_dssp_analysis_struct(dssp_list_obj, dirpath, outid)



if __name__ == '__main__':

    dirpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_nomod\rmsdmat"
    dssp_list_object_analysis_all(dirpath)