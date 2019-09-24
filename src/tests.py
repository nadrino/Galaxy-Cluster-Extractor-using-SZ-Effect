#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import healpy as hp
import matplotlib.pyplot as plt

import setup
import library

def get_file_list():
    
    files_list = []
    
    #files_list.append('LFI_SkyMap_030_1024_R2.00_full.fits')
    #files_list.append('LFI_SkyMap_044_1024_R2.00_full.fits')
    #files_list.append('LFI_SkyMap_070_2048_R2.00_full.fits')
    files_list.append('HFI_SkyMap_100_2048_R2.00_full.fits')
    files_list.append('HFI_SkyMap_143_2048_R2.00_full.fits')
    files_list.append('HFI_SkyMap_217_2048_R2.00_full.fits')
    files_list.append('HFI_SkyMap_353_2048_R2.00_full.fits')
    files_list.append('HFI_SkyMap_545_2048_R2.00_full.fits')
    files_list.append('HFI_SkyMap_857_2048_R2.00_full.fits')
    
    return files_list


def main():

    files_list = get_file_list()
    input_map, header = hp.read_map(setup.out_files_path + "Smoothed_100_GHz.fits", h=True)
    
    hp.mollview(input_map, title="Map : ", norm='hist')
    plt.show()

    print header
    
if __name__ == '__main__':
    sys.exit(main())
