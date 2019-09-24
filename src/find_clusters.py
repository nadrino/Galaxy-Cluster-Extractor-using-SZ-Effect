#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 12:32:29 2016

@author: ablanche
"""

import sys
import os
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import threading

import setup
import library

def main():

    # library.remove_dipole_and_write('ILC_SZ_comp')

    # _, background, dispersion = library.fit_pixels_distribution('ILC_SZ_comp_512', 0)
    # cleaned_map = library.clean_map('ILC_SZ_comp_512', background + 2*dispersion, 0)
    # SZ_map = hp.read_map(setup.out_files_path + 'ILC_SZ_comp_512' + ".fits")

    _, background, dispersion, d_, d_background, d_dispersion = library.fit_pixels_distribution('ILC_SZ_comp_2048', 0)
    cleaned_map = library.clean_map('ILC_SZ_comp_2048', background + 2.2*dispersion, 0)
    SZ_map = hp.read_map(setup.out_files_path + 'ILC_SZ_comp_2048' + ".fits")

    clusters = library.get_clusters_list(cleaned_map, 1e-8, SZ_map)
    print 'Nb of clusters found %s' % len(clusters)
    
    clusters.sort(key=lambda Cluster: Cluster.luminosity, reverse=True)
    library.write_clusters_info(clusters)
    
    sizes = list()
    luminosities = list()
    for i in range(len(clusters)):
        sizes.append(clusters[i].FWHM)
        luminosities.append(clusters[i].luminosity)
    
    plt.plot(sizes, luminosities, 'ro')
    plt.title('Clusters size vs luminosity')
    plt.xlabel('Size (Aperture in degrees)')
    plt.ylabel('Luminosity (in MJy / sr)')
    library.save_figure('cluster_sizes_vs_luminosity')
    
    for i in range(len(clusters)):
        cluster_name = 'Cluster_#' + str(i)
        hp.gnomview(SZ_map, rot=[clusters[i].lon, clusters[i].lat], title=cluster_name)
        library.save_figure(cluster_name, sub_folder='figures/clusters')

    return 0

if __name__ == '__main__':
    sys.exit(main())