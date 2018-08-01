#!/usr/bin/python3
"""
   Single Molecule Localization code (python port of the "pix" code).
"""

import sys
import re

#--------------------------------------------------------------------------
#   Read input parameters and calibration parameters.
#
def read_cab_para(cabfn):
    """
    Read calibration parameters from the data file, with the following format:
    <label> = <para_real> [<para_imag>] +- <para_err>
    """
    cab_para = {}
    f = open(cabfn, 'rt')
    for line in f:
        line = re.sub('\n$', '', line)
        key, line = re.split(' *= *', line)
        if (line.find('+-') > 0):
            line, err = re.split(' *\+- *', line)
        else:
            err = 0.0
        if (line.find(' ') > 0):
            v_re, v_im = re.split(' +', line);
        else:
            v_re = line
            v_im = 0.0
        var = {}
        var['name'] = key
        var['real'] = float(v_re)
        var['imag'] = float(v_im)
        var['err']  = float(err)
        cab_para[key] = var.copy()
    f.close()
    return cab_para

def read_inputfile():
    """
    Read parameters from the input file.
    """
    if (len(sys.argv) <= 1):
        print('Usage: %s: <input_file>' % sys.argv[0])
        sys.exit(0)
    inpfn = sys.argv[1]
    para  = {}

    f = open(inpfn, 'rt')
    for line in f:
        line = re.sub('[ \t]*#.*', '', line)
        line = re.sub('\n$', '', line)
        if (len(line) <= 0):
            continue

        key, val = re.split(':[ \t]*', line)
        if (key == 'IMG_TYPE'):
            para['imgfmt'] = val
        elif (key == 'INPF_RAW_IMAGE'):
            para['inpfn_img'] = val
        elif (key == 'INPF_CAB_PARA'):
            para['inpfn_cab'] = val
        elif (key == 'OUTF_SPOT'):
            para['outfn_spot'] = val
        elif (key == 'OUTF_SPOTH'):
            para['outfn_spotH'] = val
        elif (key == 'OUTF_FSTS'):
            para['outfn_fsts'] = val
        elif (key == 'OUTF_FSUM'):
            para['outfn_fsum'] = val
        elif (key == 'OUTF_CANDIDATE'):
            para['outfn_pcand'] = val
        elif (key == 'N_DIM_SML'):
            para['SML_ndim'] = val
        elif (key == 'FRAME_ID_START'):
            para['frameID1'] = int(val)
        elif (key == 'FRAME_ID_END'):
            para['frameID2'] = int(val)
        elif (key == 'FRAME_X1'):
            para['frame_x1'] = int(val)
        elif (key == 'FRAME_X2'):
            para['frame_x2'] = int(val)
        elif (key == 'FRAME_Y1'):
            para['frame_y1'] = int(val)
        elif (key == 'FRAME_Y2'):
            para['frame_y2'] = int(val)
        elif (key == 'X_FIND_PIXELS'):
            para['x_find_pixels'] = int(val)
        elif (key == 'Y_FIND_PIXELS'):
            para['y_find_pixels'] = int(val)
        elif (key == 'I_THRES_MIN'):
            para['I_thres_min'] = float(val)
        elif (key == 'I_THRES_MAX'):
            para['I_thres_max'] = float(val)
        elif (key == 'NFSEP_IDP'):
            para['nfsep'] = int(val)
        elif (key == 'NM_PIXEL_X'):
            para['nm_px_x'] = float(val)
        elif (key == 'NM_PIXEL_Y'):
            para['nm_px_y'] = float(val)
        elif (key == 'I_TO_N_PHOTONS'):
            para['i_photon'] = float(val)
        elif (key == 'RUN_MODE'):
            para['rmode'] = val
        elif (key == 'SCAN_ALGO'):
            para['scan_algo'] = int(val)
        elif (key == 'FIT_MAX_DX'):
            para['max_dx'] = float(val)
        elif (key == 'FIT_MAX_DY'):
            para['max_dy'] = float(val)
        elif (key == 'FIT_MAX_DW'):
            para['max_dwx'] = float(val)
            para['max_dwy'] = float(val)
        elif (key == 'MIN_SN_RATIO'):
            para['min_SN'] = float(val)
        elif (key == 'MIN_DI_I'):
            para['max_dI_I'] = float(val)
        elif (key == 'SOLVER_Z1'):
            para['solver_z1'] = float(val)
        elif (key == 'SOLVER_Z2'):
            para['solver_z2'] = float(val)
        elif (key == 'SOLVER_DZ'):
            para['solver_dz'] = float(val)
        elif (key == 'VERBOSE'):
            para['verb'] = int(val)
    f.close()

    if (para['imgfmt'] != 'TIFF' and para['imgfm'] != 'RAW'):
        print('!!! %s: Invalid value of parameter: IMG_TYPE' % inpfn)
        sys.exit(1)
    if (para['SML_ndim'] != '2D' and para['SML_ndim'] != '3D'):
        print('!!! %s: Invalid value of parameter: N_DIM_SML' % inpfn)
        sys.exit(1)
    if (para['rmode'] != 'SCAN' and para['rmode'] != 'SML'):
        print('!!! %s: Invalid value of parameter: RUN_MODE' % inpfn)
        sys.exit(1)
    if (para['scan_algo'] != 0 and para['scan_algo'] != 1):
        print('!!! %s: Invalid value of parameter: SCAN_ALGO' % inpfn)
        sys.exit(1)
    if (para['SML_ndim'] == '3D'):
        para['cab_para'] = read_cab_para(para['inpfn_cab'])
    return para


#--------------------------------------------------------------------------
#   Main program.
#
if __name__ == '__main__':
    para = read_inputfile()

