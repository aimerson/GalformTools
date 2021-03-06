#! /usr/bin/env python
 

import sys,os,fnmatch
import numpy as np
from .hdf5 import HDF5
from .GalformError import FileError
from .cosmology import Cosmology
from .utils.datatypes import getDataType

class GalformHDF5(HDF5):
    
    def __init__(self,*args,**kwargs):
        classname = self.__class__.__name__
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name

        # Initalise HDF5 class
        super(GalformHDF5, self).__init__(*args,**kwargs)

        # Check file is complete and not corrupted
        self.complete = True
        if "CompletionFlag" in self.fileObj.keys():
            self.complete = bool(np.array(self.fileObj["CompletionFlag"])[()])
            if not self.complete:
                raise FileError(classname+"(): File corrupted or not complete!")
        else:
            print("WARNING! "+classname+"(): cannot locate completion flag. Proceed with caution!")

        # Store version
        self.branch = np.array(self.fileObj["Version"]['branchname'])[()]
        self.version = ".".join(list(map(str,list(np.array(self.fileObj["Version"]['iversion'])))))
        self.version = self.version + "-" + np.array(self.fileObj["Version"]['revname'])[()]
        
        # Read parameters
        self.parameters = {}
        for param in self.fileObj["Parameters"].keys():
            value = np.array(self.fileObj["Parameters"][param])
            if np.ndim(value):
                self.parameters[str(param)] = np.copy(value[()])
            else:
                self.parameters[str(param)] = np.copy(value[()])
 
        # Set cosmology
        omega0 = self.parameters["omega0"]
        lambda0 = self.parameters["lambda0"]
        omegab = self.parameters["omegab"]
        h0 = self.parameters["h0"]
        sigma8 = self.parameters["sigma8"]
        ns = 1
        self.cosmology = Cosmology(omega0=omega0,lambda0=lambda0,omegab=omegab,\
                                       h0=h0,sigma8=sigma8,ns=ns)

        # Store outputs
        nout = np.array(self.fileObj["Output_Times"]["nout"])[()]
        dtype = [("a",float),("z",float),("t",float),("tlbk",float)]
        self.outputs = np.zeros(nout,dtype=dtype).view(np.recarray)
        self.outputs.a = np.array(self.fileObj["Output_Times"]["aout"])
        self.outputs.z = np.array(self.fileObj["Output_Times"]["zout"])
        self.outputs.t = np.array(self.fileObj["Output_Times"]["tout"])
        self.outputs.tlbk = np.array(self.fileObj["Output_Times"]["tlkbk"])

        return
    

    def selectOutput(self,z,returnPath=False):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        iz = np.argmin(np.fabs(self.outputs.z-z))
        outstr = "Output"+str(iz+1).zfill(3)
        if returnPath:
            out = outstr
        else:
            out = self.fileObj[outstr]
        return out
    
    def nearestRedshift(self,z):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        iz = np.argmin(np.fabs(self.outputs.z-z))
        return self.outputs.z[iz]

    def availableProperties(self,z):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        out = self.selectOutput(z,returnPath=True)
        return list(map(str,self.lsDatasets(out)))
    

    def readGalaxies(self,z,props=None):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        # Select epoch closest to specified redshift
        out = self.selectOutput(z)
        # Set list of all available properties
        allprops = self.availableProperties(z)
        # Get number of galaxies
        ngals = len(np.array(out[allprops[0]]))
        # Read all properties if not specified
        if props is None:
            props = allprops
        # Construct datatype for galaxy properties to read
        dtype = []
        for p in props:
            if len(fnmatch.filter(allprops,p))>0:
                matches = fnmatch.filter(allprops,p)
                dtype = dtype + [ (str(m),getDataType(out[m])) for m in matches ]
        galaxies = np.zeros(ngals,dtype=dtype)
        # Extract galaxy properties
        for p in galaxies.dtype.names:
            if p in allprops:                
                galaxies[p] = np.copy(np.array(out[p]))
                # Correct units for specific properties:
                # i) emission lines (correct to units of erg/s)
                if any([fnmatch.fnmatch(p,"L_"+comp+"_*") for comp in "tot total bul bulge disk".split()]):
                    galaxies[p] *= 1.0e40

        return galaxies
