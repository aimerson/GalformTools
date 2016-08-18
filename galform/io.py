#! /usr/bin/env python
 

import sys,os,fnmatch
import numpy as np


from .hdf5 import HDF5
from .GalformError import FileError


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
 
        # Store outputs
        nout = np.array(self.fileObj["Output_Times"]["nout"])[()]
        dtype = [("a",float),("z",float),("t",float),("tlbk",float)]
        self.outputs = np.zeros(nout,dtype=dtype).view(np.recarray)
        self.outputs.a = np.array(self.fileObj["Output_Times"]["aout"])
        self.outputs.z = np.array(self.fileObj["Output_Times"]["zout"])
        self.outputs.t = np.array(self.fileObj["Output_Times"]["tout"])
        self.outputs.tlbk = np.array(self.fileObj["Output_Times"]["tlkbk"])
        
        return
    

    def selectOutput(self,z):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        iz = np.argmin(np.fabs(self.outputs.z-z))
        outstr = "Output"+str(iz+1).zfill(3)
        return self.fileObj[outstr]

