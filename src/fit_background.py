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

def main():

    _, background, dispersion, error_maxvalue, error_background, error_dispersion =\
        library.fit_pixels_distribution('ILC_SZ_comp_2048', 1)
    print 'Background level %s' % background
    print 'Sigma = %s ± %s' % (dispersion, error_dispersion)

    return 0

if __name__ == '__main__':
    sys.exit(main())