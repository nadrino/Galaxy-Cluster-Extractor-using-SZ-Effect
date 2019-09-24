#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

import setup
import library

def convert_K_to_MJy(in_map, map_id):
    in_map *= setup.K_to_MJy_coef[map_id]

def convert_MJy_to_K(in_map, map_id):
    in_map /= setup.K_to_MJy_coef[map_id]


def main():

    K_map_ID = [0, 1, 2, 3]
    MJy_map_ID = [4, 5]

    # Convert K to MJy
    for i in K_map_ID:
        in_file = 'Smoothed_' + setup.frequencies[i] + 'GHz'
        print 'Converting : %s' % in_file
        in_map = hp.read_map(setup.out_files_path + in_file + '.fits')
        convert_K_to_MJy(in_map, i)
        print 'Writing converted file...'
        hp.write_map(setup.out_files_path + in_file + "_MJy" + ".fits", in_map)
        hp.mollview(in_map, title="Converted Map : " + in_file, norm='hist')
        library.save_figure(in_file, sub_folder='figures/smoothed_and_converted')

    # Already converted maps
    for i in MJy_map_ID:
        in_file = 'Smoothed_' + setup.frequencies[i] + 'GHz'
        print 'Already converted : %s' % in_file
        in_map = hp.read_map(setup.out_files_path + in_file + '.fits')
        print 'Writing file...'
        hp.write_map(setup.out_files_path + in_file + "_MJy" + ".fits", in_map)
        hp.mollview(in_map, title="Converted Map : " + in_file, norm='hist')
        library.save_figure(in_file, sub_folder='figures/smoothed_and_converted')

    plt.cla()
    plt.show()

    return 0

if __name__ == '__main__':
    sys.exit(main())