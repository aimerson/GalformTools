#! /usr/bin/env python


import sys,os,fnmatch
import numpy as np


from .hdf5 import HDF5



class GalformHDF5(HDF5):
    
    def __init__(self,*args,**kwargs):
        classname = self.__class__.__name__
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name

        # Initalise HDF5 class
        super(GalformHDF5, self).__init__(*args,**kwargs)

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
                

        return
