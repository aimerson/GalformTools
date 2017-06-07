#! /usr/bin/env python

import fnmatch,glob
import numpy as np
from ..hdf5 import HDF5
from ..cosmology import Cosmology

class LightconeHDF5(HDF5):
    
    def __init__(self,*args,**kwargs):
        classname = self.__class__.__name__
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        # Initalise HDF5 class
        super(GalformHDF5, self).__init__(*args,**kwargs)
        # Read version information
        self.galformVersionInformation = {}
        keys = self.lsDatasets("/Version/GALFORM/")
        for key in keys:
            self.galformVersionInformation[key] = np.array(self.fileObj["Version/Galform/"+key])            
        self.lightconeVersionInformation = {}
        keys = self.lsDatasets("/Version/lightcone/")
        for key in keys:
            param = np.array(self.fileObj["Version/lightcone/"+key])            
            if len(param) == 1:
                self.lightconeVersionInformation[key] = np.copy(param[0])
            else:
                self.lightconeVersionInformation[key] = np.copy(param)
        # Read parameters from header information
        self.parameters = {}
        keys = self.lsDatasets("/Header")
        for key in keys:
            self.parameters[key] = np.array(self.fileObj["Header"+key])
        # Set cosmology
        self.cosmology = Cosmology(omega0=self.parameters["Omega0"],lambda0=self.parameters["Lambda0"],\
                                       omegab=self.parameters["Omegab"],h0=self.parameters["h0"])
        # Store datatypes of each galaxy property
        self.galaxyDataTypes = {}
        if self.parameters["n_galaxies"]>0:
            for datasetName in self.lsDatasets("Data"):
                dtype = np.array(self.fileObj["Data/"+datasetName]).dtype
                self.galaxyDataTypes[datasetName] = (datasetName,dtype)
        return
    
    def availableGalaxyDatasets(self):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        return self.lsDatasets("Data")
        
    def _buildDataType(self,props):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        dummy = [self._propertiesDataTypes.append(self.galaxyDataTypes[prop]) for prop in props]
        return

    def _extractGalaxyDataset(self,datasetName):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        self._galaxyDatasets[datasetName] = np.array(self.fileObj["Data/"+datasetName])
        return

    def readGalaxies(self,props=None):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        # Construct datatype for selected properties
        if props is None:            
            self._propertiesDataType = self.galaxyDataTypes
        else:
            self._propertiesDataTypes = []
        dummy = [self._buildDataType(fnmatch.filter(self.availableGalaxyDatasets(),prop)) for prop in props]
        # Create array to store galaxy data
        self._galaxyDatasets = np.zeros(self.paramters["n_galaxies"],dtype=self._propertiesDataTypes)
        dummy = [self._extractGalaxyDataset(name) for name in self._galaxyDatasets.dtype.names]
        data = np.copy(self._galaxyDatasets)
        self._galaxyDatasets = None
        self._propertiesDataTypes = []
        return data
    

                                       

class PixelLightcone(object):
    
    def __init__(self,directory,filePrefix=None):
        classname = self.__class__.__name__
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        self.directory = directory
        # Determine prefix if not specified
        if filePrefix is None:
            filePrefix = glob.glob(self.directory+"/*.pixels-information")[0]
            filePrefix = filePrefix.replace(self.directory+"/","").replace(".pixels-information","")
        self.filePrefix = filePrefix            
        # Read in pixel information
        pixelInformationFile = self.directory+"/"+self.filePrefix+".pixels-information"
        self.pixels = np.zeros(1,dtype=[("NSIDE",int),("area",float),("number":int)]).view(np.recarray)
        with open(pixelInformationFile,'r') as f:
            line = f.readline()
            while line.startswith("#"):
                MATCH = re.search("# NSIDE = (\d+)?\n?",line)
                if MATCH is not None:
                    self.pixels.NSIDE = int(MATCH.group(1))
                MATCH = re.search("# PIXEL_AREA = ([\d\.]+) sq.deg.?\n?",line)
                if MATCH is not None:
                    self.pixels.area = float(MATCH.group(1))
                line = f.readline()
        self.pixels.number = len(np.loadtxt(pixelInformationFile,usecols=[0],dtype=int))
        # Read galaxy property locations
        propertyLocationsFile = self.directory+"/"+self.filePrefix+".property-locations"
        properties = np.loadtxt(propertyLocationsFile,dtype=[("name","|S30"),("location","|S30")],usecols=[0,1])
        self.properties = {}
        for prop,location in zip(properties["name"],properties["location"]):
            self.properties[prop] = location
        del properties
        return

    

    



