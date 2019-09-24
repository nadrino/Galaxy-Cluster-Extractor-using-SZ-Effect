#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import threading
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.stats import chisquare

import setup


class maps_processing_factory(threading.Thread):

    def __init__(self, threadID, input_file, smooth_sigma, nside_out, display=0):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.input_file = input_file
        self.smooth_sigma = smooth_sigma
        self.smooth_fwhm = self.smooth_sigma*(2*np.sqrt(2*np.log(2)))
        self.nside_out = nside_out
        
        self.input_map = hp.read_map(setup.in_files_path + input_file)
        self.smoothed_map = None
        self.smoothed_degraded_map = None

    def run(self):
        print 'Starting Thread-%s' % self.threadID
        self.smoothed_map = hp.smoothing(self.input_map, fwhm=np.radians(self.smooth_fwhm))
        self.smoothed_degraded_map = hp.ud_grade(self.smoothed_map, nside_out=self.nside_out, order_in='RING', pess=True)
        self.write_smoothed_map('Smoothed_' + setup.frequencies[self.threadID] +
                                'GHz' + '.fits')
        self.write_smoothed_and_degraded_map('Smoothed_' + setup.frequencies[self.threadID] +
                                             'GHz_' + str(self.nside_out) + 'N' + '.fits')

        print 'Thread-%s Finished' % self.threadID

    def write_smoothed_map(self, out_file_name):
        print 'Writing smoothed map in %s%s' % (setup.out_files_path, out_file_name)
        hp.write_map(setup.out_files_path + out_file_name, self.smoothed_map)

    def write_smoothed_and_degraded_map(self, out_file_name):
        print 'Writing smoothed and degraded map in %s%s' % (setup.out_files_path, out_file_name)
        hp.write_map(setup.out_files_path + out_file_name, self.smoothed_degraded_map)


def compute_S_matrix(ma_maps, F_matrix, e_vector):

    Trans_F_Matrix = np.transpose(F_matrix)
    cov = np.cov(ma_maps)
    inv_cov = np.linalg.inv(cov)

    print "\n Covariance Matrix : \n %s \n" % cov

    inv_cov = np.linalg.inv(cov)
    print np.dot(inv_cov, cov)

    print "\n Covariance Inverse Matrix : \n %s \n" % inv_cov

    S = np.dot(inv_cov, Trans_F_Matrix)
    S = np.dot(F_matrix, S)
    S = np.linalg.inv(S)
    S = np.dot(e_vector, S)
    S = np.dot(S, F_matrix)
    S = np.dot(S, inv_cov)
    print "Component : %s" % S
    S = np.dot(S, ma_maps)

    return S


def get_original_map(frequency_index, NSIDE=2048):

    #maps = []
    # truncated_maps = []
    #for i in range(len(setup.frequencies)):
        # print i
    #    maps.append(hp.read_map(setup.out_files_path + 'Smoothed_' + setup.frequencies[i] + 'GHz_MJy.fits'))

    #return maps
    str_nside = ''
    if NSIDE != 2048:
        str_nside = str(NSIDE) + 'N_'
    return hp.read_map(setup.out_files_path + 'Smoothed_' + setup.frequencies[frequency_index] +
                       'GHz_' + str_nside + 'MJy.fits')


def get_truncated_maps(NSIDE=2048):

    mask = hp.read_map(setup.out_files_path + 'mask.fits')

    ma_maps = list()
    for i in range(len(setup.frequencies)):
        ma_maps.append(np.ma.array(get_original_map(i, NSIDE), mask=(mask == 0)))

    # maps = get_original_maps()

    # mask = hp.read_map(setup.out_files_path + 'mask.fits')
    # maps = np.where(mask, maps, np.nan)
    # ma_maps = np.ma.array(maps, mask=np.tile(mask == 0, (6, 1)))
    # truncated_maps = hp.ma(maps)

    return ma_maps


def gaussian_function(x, amplitude, mean, sigma):
    """
    Compute the value of a gaussian function at x with the given parameters
    :param x:
    :param amplitude:
    :param mean:
    :param sigma:
    :return: the value of the gaussian at x
    """
    y = amplitude * np.exp(-((x - mean)**2) / (2 * sigma * sigma))
    return y


def fit_pixels_distribution(file_name, show_plot=1):

    print '-> FIT PIXELS DISTRIBUTION'

    threshold = 0.
    this_map = hp.read_map(setup.out_files_path + file_name + ".fits")
    mask = hp.read_map(setup.out_files_path + "mask.fits")

    S_mask = np.ma.array(this_map, mask=(mask == 0), fill_value=hp.UNSEEN)

    bin_number = 8000
    array_map = np.array(S_mask)
    y, x = np.histogram(array_map, bin_number, range=(-5e-5, 1e-4))
    max_y = np.float(np.max(y))
    normal_y = y/max_y
    max_x = np.float(np.max(x))
    normal_x = x[:-1]/max_x
    fit_parameters, covariant = curve_fit(gaussian_function, normal_x, normal_y)
    # chisquare(y, ...)

    print 'Best fit = %s' % fit_parameters
    error = np.sqrt(np.diag(covariant))
    print 'Error = %s' % error

    print 'Covariance matrix :'
    print covariant

    y = np.append(y,0.)

    maxvalue = fit_parameters[0] * max_y
    error_maxvalue = error[0] * max_y
    background = fit_parameters[1] * max_x
    error_background = error[1] * max_x
    dispersion = np.abs(fit_parameters[2] * max_x)
    error_dispersion = np.abs(error[2] * max_x)
    print 'Max = %s, Mean = %s, Dispersion = %s' % (maxvalue, background, dispersion)

    if show_plot:

        plt.plot(x, y, 'b+:', label='data')
        plt.plot(x, gaussian_function(x, maxvalue, background, dispersion),
                     'r.:', label='fit')
        plt.legend()
        plt.title('Flux distribution')
        plt.xlabel('Amplitude')
        plt.ylabel('Frequency of appearance')

        plt.show()

    return maxvalue, background, dispersion, error_maxvalue, error_background, error_dispersion


def remove_dipole_and_write(file_name):

    print '-> REMOVING DIPOLE'
    ILC_SZ_map = hp.read_map(setup.out_files_path + file_name + '.fits')
    hp.mollview(ILC_SZ_map, title='BEFORE REMOVING DIPOLE', norm='hist')
    hp.remove_dipole(ILC_SZ_map)
    hp.mollview(ILC_SZ_map, title='AFTER REMOVING DIPOLE', norm='hist')
    hp.write_map(setup.out_files_path + file_name + '_dipole_rem', ILC_SZ_map)
    plt.show()

    return 0


def clean_map(file_name, threshold, show_plot=1):

    this_map = hp.read_map(setup.out_files_path + file_name + ".fits")
    mask = hp.read_map(setup.out_files_path + "mask.fits")
    cleaned_map = np.ma.array(this_map, mask=(mask == 0), fill_value=0.)
    cleaned_map = (cleaned_map - threshold)*(cleaned_map > threshold)#*(cleaned_map < 1e19)
    if show_plot:
        hp.mollzoom(cleaned_map, title="Cleaned Map : " + file_name)#, min = 0, max = 1e-4)
        plt.show()

    return cleaned_map


class Cluster:

    def __init__(self, pixel_index, pixel_value):

        self.pixel_list = []
        self.pixel_list.append(pixel_index)

        self.peak_index = pixel_index
        self.peak_value = pixel_value
        self.integral = pixel_value

        self.lon = 0.
        self.lat = 0.

        self.pixel_angle = 0.
        self.FWHM = 0.
        self.luminosity = 0
        self.luminosity_error = 0
        self.touched_mask = 0

    def add_pixel(self, pixel_index, pixel_value):

        self.pixel_list.append(pixel_index)
        self.integral += pixel_value
        # print 'Nb %s : %s' % (len(self.pixel_list), pixel_value)
        # print np.degrees(hp.pix2ang(512, pixel_index))

        if self.peak_value < pixel_value:
            self.peak_value = pixel_value
            self.peak_index = pixel_index

    def check_neighbourhood(self, pixels_map, threshold, check_map, pixel_index):

        neighbours_indexes = hp.get_all_neighbours(hp.get_nside(pixels_map), pixel_index)
        hp.get_all_neighbours(hp.get_nside(pixels_map), neighbours_indexes[0])
        # Loop over all neighbours
        for i in range(len(neighbours_indexes)):
            # Checking if its already checked
            if check_map[neighbours_indexes[i]] == 0:
                check_map[neighbours_indexes[i]] = 1
                # Checking if the selected pixel is in the cluster
                if pixels_map[neighbours_indexes[i]] >= threshold and self.touched_mask == 0:
                    self.add_pixel(neighbours_indexes[i], pixels_map[neighbours_indexes[i]])
                    self.check_neighbourhood(pixels_map, threshold, check_map, neighbours_indexes[i])
            elif check_map[neighbours_indexes[i]] == 2:
                self.touched_mask = 1

    def compute_luminosity(self, SZ_map):

        self.compute_FWHM_in_degrees(SZ_map)
        NSIDE = hp.get_nside(SZ_map)

        pixels_in_disk = hp.query_disc(NSIDE, hp.pix2vec(NSIDE, self.peak_index), np.radians(self.FWHM*1.5))
        nb_pixels_disk = len(pixels_in_disk)
        vania_luminosity = 0

        for pixel in pixels_in_disk:
            vania_luminosity += SZ_map[pixel]

        inner_ring = hp.query_disc(NSIDE, hp.pix2vec(NSIDE, self.peak_index), np.radians(self.FWHM*3.5))
        outer_ring = hp.query_disc(NSIDE, hp.pix2vec(NSIDE, self.peak_index), np.radians(self.FWHM*5.5))
        pixels_values_ring = list()
        nb_pixels_background = len(outer_ring) - len(inner_ring)

        for pixel in outer_ring:
            if not(pixel in inner_ring):
                pixels_values_ring.append(SZ_map[pixel])

        self.luminosity_error = np.std(pixels_values_ring)
        local_noise_ring = sum(pixels_values_ring)
        self.luminosity = vania_luminosity - local_noise_ring * (nb_pixels_disk / nb_pixels_background)
        print 'luminosity = %s \pm %s' %(self.luminosity, self.luminosity_error)

    def compute_pixel_angle_in_degrees(self, NSIDE):

        self.pixel_angle = (1.7 * (2048/NSIDE)) / 60.

    def compute_FWHM_in_degrees(self, SZ_map):

        NSIDE = hp.get_nside(SZ_map)

        self.compute_pixel_angle_in_degrees(NSIDE)
        intensity_profile = list()
        intensity_profile.append(SZ_map[self.peak_index])
        half_max = SZ_map[self.peak_index]/2.
        upper_bound_aperture = 0.
        aperture = list()
        aperture.append(0.)

        '''
        Arbitrary go to to maximum disk of 1 deg (assuming no cluster is bigger)
        '''
        nb_rings = int(2. / self.pixel_angle)
        for i in range(nb_rings):
            aperture.append((i+2)*self.pixel_angle)
            pixels_in_disk = hp.query_disc(NSIDE, hp.pix2vec(NSIDE, self.peak_index), np.radians(aperture[i+1]))
            #disk_area = len(pixels_in_disk)
            ring_area = 0            
            ring_luminosity = 0.
            for pixel in pixels_in_disk:
                ring_luminosity += SZ_map[pixel]
                ring_area += 1
            new_intensity_profile = ring_luminosity/ring_area
            intensity_profile.append(new_intensity_profile)

        intensity_profile -= intensity_profile[len(intensity_profile)-1]
        upper_bound_aperture_index = 0
        for i in range(nb_rings):
            if intensity_profile[i] > half_max:# and trigger==False:
                upper_bound_aperture_index = i
            else:
                break

        # Interpolation using tales theorem
        delta_y1 = intensity_profile[upper_bound_aperture_index] - half_max
        delta_y = intensity_profile[upper_bound_aperture_index] - intensity_profile[upper_bound_aperture_index + 1]
        delta_x = aperture[upper_bound_aperture_index + 1] - aperture[upper_bound_aperture_index]
        delta_x1 = (delta_y1/delta_y)*delta_x

        self.FWHM = aperture[upper_bound_aperture_index] + delta_x1
        print 'self.FWHM %s' % self.FWHM
        # print 'Aperture FWHM = %s' % np.interp(SZ_map[self.peak_index]/2., intensity_profile, aperture)
        #plt.plot(aperture, intensity_profile)
        #plt.show()

    def compute_lon_lat(self, NSIDE):

        theta, phi = np.degrees(hp.pix2ang(NSIDE, self.peak_index))
        self.lat = 90. - theta
        if phi >= 180:
            self.lon = -(360 - phi)
        else:
            self.lon = phi


def get_clusters_list(pixels_map, threshold=0., SZ_map=None):

    print "-> GETTING CLUSTERS"

    clusters_list = []
    size = pixels_map.shape[0]
    mask = hp.read_map(setup.out_files_path + "mask.fits")
    check_map = (mask < 1)*2

    for pixel_index in range(size):
        if check_map[pixel_index] == 0:
            check_map[pixel_index] = 1
            if pixels_map[pixel_index] > threshold:
                new_cluster = create_new_cluster(pixels_map, threshold, check_map, pixel_index)
                if new_cluster.touched_mask == 0 and len(new_cluster.pixel_list) > 30:
                    print 'CLUSTER #%s FOUND !' % len(clusters_list)
                    if SZ_map != None:
                        new_cluster.compute_luminosity(SZ_map)
                        new_cluster.compute_lon_lat(hp.get_nside(SZ_map))

                    clusters_list.append(new_cluster)

    print "-> END GETTING CLUSTERS"
    return clusters_list


def create_new_cluster(pixels_map, threshold, check_map, first_pixel_index):

    # Recursion security
    sys.setrecursionlimit(10000)
    new_cluster = Cluster(first_pixel_index, pixels_map[first_pixel_index])

    new_cluster.check_neighbourhood(pixels_map, threshold, check_map, first_pixel_index)

    return new_cluster


def theta_phi_to_lon_lat(theta, phi):

    lat = 90 - theta
    if phi >=180:
        lon = -(360 - phi)
    else:
        lon = phi

    return lon, lat


def save_figure(out_file_name, sub_folder='figures'):
    
    try:
        os.mkdir(setup.out_files_path + 'figures/')
    except OSError:
        _ = 0

    folder_path = setup.out_files_path + sub_folder + '/'
    try:
        os.mkdir(folder_path)
    except OSError:
        _ = 0

    print 'Writing figure...'
    plt.savefig(folder_path + out_file_name + '.png')
    print 'PNG file "%s" saved in %s' % (out_file_name + '.png', folder_path)


def write_clusters_info(clusters_list):

    txt_file = open(setup.out_files_path + "clusters.txt", "w")
    for i in range(len(clusters_list)):
        string = str(clusters_list[i].lon) + ' ' +\
                 str(clusters_list[i].lat) + ' ' +\
                 str(clusters_list[i].FWHM) + ' ' +\
                 str(clusters_list[i].luminosity) + ' ' +\
                 str(clusters_list[i].luminosity_error) + '\n'
        txt_file.write(string)

    txt_file.close()
