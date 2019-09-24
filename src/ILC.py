#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 14:20:45 2016

@author: ablanche
"""

import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

import setup
import library

Total_matrix = [setup.SZ_comp, setup.K_to_MJy_coef, setup.dust_comp, setup.CO_comp, setup.synchrotron_comp]

def main():

    nb_comp = 2
    NSIDE = 2048
    F_matrix = []

    for i in range(nb_comp):
        F_matrix.append(Total_matrix[i])

    print 'F_matrix : %s' % F_matrix

    truncated_maps = library.get_truncated_maps(NSIDE)
    mask = hp.read_map(setup.out_files_path + 'mask.fits')

    for i in range(nb_comp):

        #if i == 0:
        if True:

            print 'Computing component %s' % i

            e_vector = list()

            for j in range(nb_comp):
                e_vector.append(0)

            e_vector[i] = 1.

            S = library.compute_S_matrix(truncated_maps, F_matrix, e_vector)
            # hp.mollview(S[i], title="Component %s" % i, norm='hist')

            if i == 0:
                hp.gnomview(S, rot=[0, 90])
                hp.mollzoom(S, norm='hist')
                hp.write_map(setup.out_files_path + "ILC_SZ_comp_" + str(NSIDE) + ".fits", S)

            S = np.ma.array(S, mask=(mask==0), fill_value=hp.UNSEEN)
            hp.mollview(S, title="Component_%s" % i, norm='hist')
            library.save_figure('Component_%s' % i, sub_folder='figures/ILC')
            del S

    plt.show()

    return 0

if __name__ == '__main__':
    sys.exit(main())
