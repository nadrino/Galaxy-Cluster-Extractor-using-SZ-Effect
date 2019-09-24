#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np

import setup
import library

def get_file_list():

    files_list = []

    #files_list.append('LFI_SkyMap_030_1024_R2.00_full.fits')
    #files_list.append('LFI_SkyMap_044_1024_R2.00_full.fits')
    #files_list.append('LFI_SkyMap_070_2048_R2.00_full.fits')
    files_list.append('HFI_SkyMap_100_2048_R2.02_full.fits')
    files_list.append('HFI_SkyMap_143_2048_R2.02_full.fits')
    files_list.append('HFI_SkyMap_217_2048_R2.02_full.fits')
    files_list.append('HFI_SkyMap_353_2048_R2.02_full.fits')
    files_list.append('HFI_SkyMap_545_2048_R2.02_full.fits')
    files_list.append('HFI_SkyMap_857_2048_R2.02_full.fits')

    return files_list


def main():

    files_list = get_file_list()

    # Getting all PSF sigmas in arcmin]
    frequencies = ['100', '143', '217', '353', '545', '857']
    #FWHMs = [32.3, 27, 13.2, 9.7, 7.3, 5, 4.9, 4.8, 4.6]
    FWHMs = [9.7, 7.3, 5, 4.9, 4.8, 4.6]
    sigmas = FWHMs / (2*np.sqrt(2*np.log(2)))

    # Convert [arcmin] to [degrees]
    sigmas /= 60.

    # Computing maps in parallel
    threads = []
    for i in range(len(files_list)):
        smooth_sigma = np.sqrt(sigmas[0]**2 - sigmas[i]**2)
        threads.append(library.maps_processing_factory(i, files_list[i], smooth_sigma, 512))
        threads[i].start()



    return 0


if __name__ == '__main__':
    sys.exit(main())
