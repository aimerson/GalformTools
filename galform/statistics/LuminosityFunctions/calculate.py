#! /usr/bin/env python


import sys
import copy
import numpy as np
from ...utils.progress import Progress
from ...io import GalformHDF5



class ComputeLuminosityFunction(object):
    
    def __init__(self,magnitudeBins=None,luminosityBins=None):        
        classname = self.__class__.__name__
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        # Store bins for magnitudes/luminosities
        self.magnitudeBins = magnitudeBins
        if self.magnitudeBins is None:
            self.magnitudeBins = np.arange(-40.0,-5.0,0.2)
        self.luminosityBins = luminosityBins
        if self.luminosityBins is None:
            self.luminosityBins = np.linspace(35.0,44.0,100)
        # Dictionary to store results
        self.luminosityFunction = {}
        return

    
    def processOutput(self,galaxiesFile,redshifts=None,props=None):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        # Open GALFORM file
        GALFORM = GalformHDF5(galaxiesFile,'r')
        # Set redshifts to process
        if redshifts is None:
            redshifts = np.copy(GALFORM.outputs.z)
        # Loop over redshifts
        for z in redshifts:            
            # Select properties
            allProperties = GALFORM.availableProperties(z)
            if props is None:
                goodProps = allProperties
            else:
                goodProps = []
                for p in tmpProps:
                    goodProps = goodProps + fnmatch.filter(allProps,p)
            # Select output
            outstr = selectOutput(z,returnPath=False)
            out = selectOutput(z,returnPath=True)
            # Loop over properties
            redshiftLF = {}
            PROG = Progress(len(goodProps))
            for p in goodProps:
                values = np.array(out[p])
                if fnmatch.fnmatch(p,"L_*") or fnmatch.fnmatch(p,"Lmod_*") or fnmatch.fnmatch(p,"Ld_*[or]"):
                    values np.log10(values) + 40.0
                    bins = self.luminosityBins
                elif fnmatch.fnmatch(p,"Ld_*"):
                    bins = self.magnitudeBins
                else:
                    bins = None
                    values = None
                redshiftLF[p],bins = np.histogram(values,bins=bins,weights=weight)
                PROG.increment()
                if verbose:
                    PROG.print_status_line()
            self.luminosityFunction[outstr] = copy.copy(redshiftLF)
        GALFORM.close()
        return 

    
    def addLuminosityFunctions(self,lfObj):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        if verbose:
            print(funcname+"(): adding luminosity functions...")
        # Check whether luminosity functions for this class have
        # already been calculated.  If not simply store second object
        # as luminosity functions for this class.
        if len(self.luminosityFunction.keys()) == 0:
            self.luminosityFunction = lfClass.luminosityFunction.copy()
            self.magnitudeBins = np.copy(lfClass.magnitudeBins)
            self.luminosityBins = np.copy(lfClass.luminosityBins)
        else:
            # Check luminosity and magnitude bins for two LF objects
            # are consistent
            binsDiff = np.fabs(self.luminosityBins-lfObj.luminosityBins)
            if any(binsDiff>binsTolerance):
                raise ValueError(funcname+"(): cannot add luminosity functions -- luminosity bins are not consistent!")
            binsDiff = np.fabs(self.magnitudeBins-lfObj.magnitudeBins)
            if any(binsDiff>binsTolerance):
                raise ValueError(funcname+"(): cannot add luminosity functions -- magnitude bins are not consistent!")
            # Add luminosity functions
            PROG = Progress(len(self.luminosityFunction.keys()))
            for outKey in self.luminosityFunction.keys():
                if outKey in lfObj.luminosityFunction.keys():
                    for p in self.luminosityFunction[outKey].keys():
                        if p in lfObj.luminosityFunction[outKey].keys():
                            self.luminosityFunction[outKey][p] += lfClass.luminosityFunction[outKey][p]
                PROG.increment()
                if verbose:
                    PROG.print_status_line()
        return


      def writeToHDF5(self,hdf5File,verbose=False):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        if verbose:
            print(funcname+"(): Writing luminosity function data to "+hdf5File+" ...")
        fileObj = HDF5(hdf5File,'w')
        # Write luminosity and magnitude bins
        luminosityBinWidth = self.luminosityBins[1] - self.luminosityBins[0]
        luminosityBins = self.luminosityBins[:-1] + luminosityBinWidth/2.0
        fileObj.addDataset("/","luminosityBins",luminosityBins,chunks=True,compression="gzip",\
                               compression_opts=6)
        fileObj.addAttributes("/luminosityBins",{"units":"log10(erg/s)"})
        magnitudeBinWidth = self.magnitudeBins[1] - self.magnitudeBins[0]
        magnitudeBins = self.magnitudeBins[:-1] + magnitudeBinWidth/2.0
        fileObj.addDataset("/","magnitudeBins",magnitudeBins,chunks=True,compression="gzip",\
                               compression_opts=6)
       # Create outputs group and write data for each output in luminosity functions dictionary                                                                                                                            
        fileObj.mkGroup("Outputs")
        for outstr in self.luminosityFunction.keys():
            fileObj.mkGroup("Outputs/"+outstr)
            iout = int(outstr.replace("Output",""))
            iselect = np.argwhere(self.galHDF5Obj.outputs.iout==iout)
            z = str(self.galHDF5Obj.outputs["z"][iselect][0][0])
            fileObj.addAttributes("Outputs/"+outstr,{"redshift":z})
            for p in self.luminosityFunction[outstr].keys():
                path = "Outputs/"+outstr+"/"
                redshiftLabel = fnmatch.filter(p.split(":"),"z*")[0]
                lfData = self.luminosityFunction[outstr][p]
                if fnmatch.fnmatch(p,"*LineLuminosity*"):
                    lfData /= luminosityBinWidth
                else:
                    lfData /= magnitudeBinWidth
                fileObj.addDataset(path,p.replace(":"+redshiftLabel,""),lfData,chunks=True,compression="gzip",\
                                       compression_opts=6)
        fileObj.close()
        if verbose:
            print(funcname+"(): luminosity function data successfully written to "+hdf5File)
        return



class GalformLuminosityFunction(object):

    def __init__(self,luminosityFunctionFile):
        classname = self.__class__.__name__
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        self.file = luminosityFunctionFile
        f = HDF5(self.file,'r')
        # Read bins arrays
        bins = f.readDatasets("/",required=["luminosityBins","magnitudeBins"])
        self.luminosityBins = np.copy(bins["luminosityBins"])
        self.magnitudeBins = np.copy(bins["magnitudeBins"])
        del bins
        # Read list of available outputs and datasets
        self.outputs = list(map(str,f.lsGroups("Outputs")))
        self.redshifts = np.ones(len(self.outputs))
        self.datasets = {}
        for i,out in enumerate(self.outputs):
            self.redshifts[i] = f.readAttributes("Outputs/"+out)["redshift"]
            self.datasets[out] = list(map(str,f.lsDatasets("Outputs/"+out)))
        f.close()
        return

    def getDatasets(self,z,required=None,verbose=False):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        iselect = np.argmin(np.fabs(self.redshifts-z))
        path = "Outputs/"+self.outputs[iselect]
        if verbose:
            print_str = funcname+"(): Reading luminosity function(s) for z = "+\
                str(self.redshifts[iselect])+"\n         -- located in path "+path
            print(print_str)
        f = HDF5(self.file,'r')
        availableDatasets = list(map(str,f.lsDatasets(path)))
        datasets = []
        for req in required:
            datasets = datasets + fnmatch.filter(availableDatasets,req)
        lfData = f.readDatasets(path,required=datasets)
        f.close()
        return lfData

