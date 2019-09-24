#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 12:32:29 2016

@author: ablanche
"""

import sys
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import threading

import setup
import library

def create_mask(in_map, threshold):

    mask = in_map < threshold
    return mask


def main():

    in_file = 'Smoothed_857GHz_MJy'
    #in_file = 'Smoothed_857GHz_512N_MJy'
    in_map = hp.read_map(setup.out_files_path + in_file + '.fits')
    mask = create_mask(in_map, 11)
    unmasked_area = mask.sum()
    #for i in range(len(mask)):
    #    if mask[i] == 1:
    #        masked_area += 1.
    
    hp.write_map(setup.out_files_path + "mask.fits", mask)
    fig_title = "Keeping %s pourcent of the sky" % ( (unmasked_area*100.)/len(mask) )
    hp.mollview(mask, title=fig_title, norm='hist')
    library.save_figure('current_mask', sub_folder='figures/mask')

    plt.show()

    return 0

if __name__ == '__main__':
    sys.exit(main())