#! /usr/bin/env python

import copy
import sys,os,fnmatch
import numpy as np
from .io import GalformHDF5
from .utils.datatypes import getDataType
from .cosmology import Cosmology
from .constants import megaParsec,centi,Pi


class GalformEmissionLines(object):
    
    def __init__(self,galformObj):
        classname = self.__class__.__name__
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        # Set up dictionary of available emission lines
        linesObj = galformObj.fileObj['Lines']
        nlines = int(np.array(linesObj['nline'])[()])            
        dtype = []
        lineKeys = linesObj.keys()
        dummy = lineKeys.pop(lineKeys.index("nline"))
        dummy = [dtype.append((str(k),getDataType(np.array(linesObj[k])))) for k in lineKeys]
        self.emissionLines = np.zeros(nlines,dtype=dtype)
        for name in self.emissionLines.dtype.names:
            self.emissionLines[name] = np.copy(np.array(linesObj[name]))
        self.emissionLines = self.emissionLines.view(np.recarray)
        # Create instance of cosmology object
        self.cosmology = copy.copy(galformObj.cosmology)
        return

    def getWavelength(self,linename):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        index = np.argwhere(self.emissionLines.linename==linename)
        return self.emissionLines.lambda_line[index][0][0]
        

    def flux(self,luminosity,redshift):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        if np.ndim(luminosity) > 0 and type(redshift) == float:
            redshift = np.ones_like(luminosity)*redshift
        dL = self.cosmology.luminosity_distance(redshift)*megaParsec/centi
        area = 4.0*Pi*dL**2
        return luminosity/area
    
    
            
        
