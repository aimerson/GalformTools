#! /usr/bin/env python

import copy
import sys,os,fnmatch
import numpy as np
from .cosmology import Cosmology


def apparentMagnitude(absMag,redshift,cosmologyObj):
    if np.ndim(redshift) == 0:
        redshift = np.ones_like(absMag)*redshift
    appMag = absMag + cosmologyObj.band_corrected_distance_modulus(redshift)
    return appMag

def absoluteMagnitude(appMag,redshift,cosmologyObj):
    if np.ndim(redshift) == 0:
        redshift = np.ones_like(appMag)*redshift
    absMag = appMag - cosmologyObj.band_corrected_distance_modulus(redshift)
    return absMag
    

