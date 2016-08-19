#! /usr/bin/env python
 

import sys,os,fnmatch
import numpy as np


from .hdf5 import HDF5
from .GalformError import FileError
from .cosmology import Cosmology

class GalformHDF5(HDF5):
    
    def __init__(self,*args,**kwargs):
        classname = self.__class__.__name__
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name

        # Initalise HDF5 class
        super(GalformHDF5, self).__init__(*args,**kwargs)

        # Check file is complete and not corrupted
        complete = bool(np.array(self.fileObj["CompletionFlag"])[()])
        if not complete:
            raise FileError(classname+"(): File corrupted or not complete!")

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

        
        #print self.fileObj.keys()

               

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


    def availableProperties(self,z):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        out = self.selectOutput(z,returnPath=True)
        return list(map(str,self.lsDatasets(out)))
