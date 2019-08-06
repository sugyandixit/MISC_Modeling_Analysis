import os


class DSSPObject(object):

    """
    create a dssp object to store vital information from the .dssp file
    """

    def __init__(self):
        """
        initialize dssp object
        """
        self.dssp_fpath = None
        self.accessible_surface = None
        self.total_num_hbonds = None
        self.total_num_hbonds_parallel_bridge = None
        self.total_num_hbonds_antiparallel_bridge = None
        self.total_num_hbonds_i5_minus = None
        self.total_num_hbonds_i4_minus = None
        self.total_num_hbonds_i3_minus = None
        self.total_num_hbonds_i2_minus = None
        self.total_num_hbonds_i1_minus = None
        self.total_num_hbonds_i0_plus = None
        self.total_num_hbonds_i1_plus = None
        self.total_num_hbonds_i2_plus = None
        self.total_num_hbonds_i3_plus = None
        self.total_num_hbonds_i4_plus = None
        self.total_num_hbonds_i5_plus = None

        self.res_num = []
        self.aa_code = []
        self.struct_code = []
        self.bridge_partner_1 = []
        self.bridge_partner_2 = []
        self.bridge_label = []
        self.acc = []
        self.tco = []
        self.kappa = []
        self.alpha = []
        self.phi = []
        self.psi = []
        self.x_ca = []
        self.y_ca = []
        self.z_ca = []




def create_dssp_obj(dssp_fpath):
    """
    given the dssp fpath, store the information in an object
    :param dssp_fpath: dssp file path
    :return: dssp object
    """

    dssp_obj = DSSPObject()

    dssp_obj.dssp_fpath = dssp_fpath


    with open(dssp_fpath, 'r') as dsspfile:
        dssp_file = dsspfile.read().splitlines()

        for ind, line in enumerate(dssp_file):
            chars = line.split()
            if ind == 7:
                dssp_obj.accessible_surface = float(chars[0])
            if ind == 8:
                dssp_obj.total_num_hbonds = float(chars[0])
            if ind == 9:
                dssp_obj.total_num_hbonds_parallel_bridge = float(chars[0])
            if ind == 10:
                dssp_obj.total_num_hbonds_antiparallel_bridge = float(chars[0])
            if ind == 11:
                dssp_obj.total_num_hbonds_i5_minus = float(chars[0])
            if ind == 12:
                dssp_obj.total_num_hbonds_i4_minus = float(chars[0])
            if ind == 13:
                dssp_obj.total_num_hbonds_i3_minus = float(chars[0])
            if ind == 14:
                dssp_obj.total_num_hbonds_i2_minus = float(chars[0])
            if ind == 15:
                dssp_obj.total_num_hbonds_i1_minus = float(chars[0])
            if ind == 16:
                dssp_obj.total_num_hbonds_i0_plus = float(chars[0])
            if ind == 17:
                dssp_obj.total_num_hbonds_i1_plus = float(chars[0])
            if ind == 18:
                dssp_obj.total_num_hbonds_i2_plus = float(chars[0])
            if ind == 19:
                dssp_obj.total_num_hbonds_i3_plus = float(chars[0])
            if ind == 20:
                dssp_obj.total_num_hbonds_i4_plus = float(chars[0])
            if ind == 21:
                dssp_obj.total_num_hbonds_i5_plus = float(chars[0])
            if ind > 27:

                dssp_obj.res_num.append(int(line[5:10].strip()))
                dssp_obj.aa_code.append(line[13])
                dssp_obj.struct_code.append(line[16])
                dssp_obj.bridge_partner_1.append(int(line[26:29].strip()))
                dssp_obj.bridge_partner_2.append(int(line[30:33].strip()))
                dssp_obj.bridge_label.append(line[33])
                dssp_obj.acc.append(float(line[34:38].strip()))
                dssp_obj.tco.append(float(line[85:91].strip()))
                dssp_obj.kappa.append(float(line[91:97].strip()))
                dssp_obj.alpha.append(float(line[97:103].strip()))
                dssp_obj.phi.append(float(line[103:109].strip()))
                dssp_obj.psi.append(float(line[109:115].strip()))
                dssp_obj.x_ca.append(float(line[115:122].strip()))
                dssp_obj.y_ca.append(float(line[122:129].strip()))
                dssp_obj.z_ca.append(float(line[129:136].strip()))


    return dssp_obj



if __name__ == '__main__':

    dssp_fpath = r"C:\Users\sugyan\Documents\MembraneMDfiles\Serf_tmp_model\Replica_Exchange_Analysis_mod\rmsd_mat\repid_0_frame_12196.dssp"
    dssp_obj = create_dssp_obj(dssp_fpath)
    print('heho')
