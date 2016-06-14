#! /usr/bin/env python

import sys,os,fnmatch
import numpy as np
from .progress import *
from .plotting.utils import *

def effective_wavelength(wav,resp,alpha=1.0,beta=0.0):
    from scipy import integrate
    wav_eff = 0.0
    flux1 = wav**(-1.0*alpha - beta - 1.0)*resp
    flux2 = wav**(-1.0*beta - 1.0)*resp
    area2 = np.sum(integrate.cumtrapz(flux2,wav))
    if area2 > 0.0:
        area1 = np.sum(integrate.cumtrapz(flux1,wav))
        wav_eff  = (area1/area2)**(-1.0/alpha)
    return wav_eff


class GalformFilters(object):
    """
    CLASS GalformFilters

    Python class to query a GALFORM filters file.
    
    USAGE: filters = GalformFilters(filter_file)
    
    FUNCTIONS:

      ids = filters.names() -- returns list of short IDs of filters
      n = filters.number_of_filters() -- returns total number of filters
                                         in filters file
      i = filters.filter_number(id) -- given the short ID, returns the
                                       integer number of a filter
      wav,trans = filters.filter(id) -- returns wavelenght and transmission
                                        for filter with specified short ID
      ax = filters.plot(id,ax,[ls],[colour],[label]) -- plots the filter with specified short ID to
                                                        specified axis. Options: [ls] = line-style,
                                                        [colour] = line colour, [label] = label for
                                                        legend of plot.                                        
    """
    def __init__(self,fileName,**kwargs):
        classname = self.__class__.__name__
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        # Set verbosity
        if "verbose" in kwargs.keys():
            self._verbose = kwargs["verbose"]
        else:
            self._verbose = False
        # Open filters file
        self.fileName = fileName
        if self._verbose:
            print(classname+"(): FILTERS FILE = "+self.fileName)
            
        # Read in filters file
        with open(self.fileName,'r') as f:
            filecontents = f.read().split('\n')
        # Store version of filters file
        self.version = fnmatch.filter(filecontents,"*Filter response file - version:*")[0].split()[-1]
        if self._verbose:
            print(classname+"(): VERSION = "+self.version)
        # Get number of filters
        self.Nfil = int(fnmatch.filter(filecontents,"*Number of filters in file*")[0].split()[-1])
        if self._verbose:
            print(classname+"(): NUMBER OF FILTERS = "+str(self.Nfil))
        # Extract information for individual filters                        
        lines = fnmatch.filter(filecontents,"#     [0-9]*")+\
            fnmatch.filter(filecontents,"#    [0-9][0-9]*")+\
            fnmatch.filter(filecontents,"#   [0-9][0-9][0-9]*")
        # Check number of filters matches number at top of file
        if len(lines) != self.Nfil:
            print("WARNING! "+classname+"(): Number of filters in file does not match value in header!")                
            print("         Correcting number of filters. Number of filters = "+str(len(lines)))
            self.Nfil = len(lines)
            
        # Store list of filters
        dtype = [("index",int),("nbins",int),("hdrStart",int),('dataStart',int),("ID",'|S2'),("longID","|S50")]
        self.filters = np.zeros(self.Nfil,dtype=dtype).view(np.recarray)            
        if self._verbose:
            print(classname+"(): Reading filters...")
        progress = Progress(self.Nfil)
        for i in range(self.Nfil):
            l = lines[i].replace("#","")                
            self.filters["index"][i] = int(l.split()[0])
            self.filters["nbins"][i] = int(l.split()[1])
            self.filters["ID"][i] = l.split()[2]
            self.filters["longID"][i] = " ".join(l.split()[3:])            
            istart = filecontents.index(lines[i])            
            self.filters["hdrStart"][i] = istart
            while filecontents[istart].startswith("#"):
                istart += 1
            self.filters["dataStart"][i] = istart
            progress.increment()
            if self._verbose:
                progress.print_status_line()
        self.filters = self.filters.view(np.recarray)
        del filecontents
        return        

    def outputIDs(self,outfile=None):        
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        # Output list of filter IDs for reference
        if outfile is None:
            outfile = self.fileName.replace(".dat","_ids.dat")            
        f = open(outfile,"w")
        f.write("# VERSION= "+self.version+"\n")
        f.write("# NFILTERS= "+str(self.Nfil)+"\n")
        f.write("# Columns: index | nbins | header start | data start | ID | longID\n")
        fmt = "%4i | %4i | %10i | %10i | %2s | %50s"
        np.savetxt(f,self.filters,fmt=fmt,delimiter='|')
        f.close()
        if self._verbose:
            print(funcname+"(): written filter IDs to file "+outfile)
        return

    def getNumber(self,shortID):
        funcname = self.__class__.__name__+"."+sys._getframe().f_code.co_name
        iargs = np.argwhere(self.filters.ID==shortID)
        if len(iargs) == 0:
            if self._verbose:
                print("WARNING! "+funcname+"(): filter with ID '"+shortID+"' not found in file!")
            return None
        else:
            if len(iargs) == 1:
                return self.filters.index[iargs[0]][0]
            else:
                print("WARNING! "+funcname+"(): found multiple filters with ID '"+shortID+"'!")
                print("         Available filters include:")
                for iarg in iargs:
                    line = '          ( {0:3.0f} )  {1:s}'.format(self.filters.index[iarg[0]],\
                                                                      self.filters.longID[iarg[0]])
                    print(line)
                return iargs

    def getID(self,NUM,short=True):
        iarg = self.filters.index==NUM
        if short:
            return self.filters.ID[iarg][0]
        else:
            return self.filters.longID[iarg][0]
        
    def get(self,filterNum,output=None):
        import itertools
        import tempfile
        fout = None
        if output is not None:
            fout = open(output,'w')
        if not isinstance(filterNum,list):
            filterNum = [filterNum]
        progress = Progress(len(filterNum))        
        dtype = [("wavelength",float),("transmission",float)]
        filterInfo = {}
        with open(self.fileName,'r') as f:
            for NUM in filterNum:
                hstart = self.filters.hdrStart[NUM-1]                
                dstart = self.filters.dataStart[NUM-1]                    
                ndata = self.filters.nbins[NUM-1]    
                header = "".join(list(itertools.islice(f,hstart,dstart)))
                if fout is not None:
                    fout.write(header)
                f.seek(0)
                data = "".join(list(itertools.islice(f,dstart,dstart+ndata)))
                if fout is not None:
                    fout.write(data)
                with tempfile.NamedTemporaryFile() as temp:
                    temp.write(data)
                    temp.flush()                
                    temp.seek(0)
                    data = np.loadtxt(temp,dtype=dtype)
                filterInfo[str(NUM)] = (header,np.copy(data).view(np.recarray))
        if fout is not None:
            fout.close()
        return filterInfo


    def plot(self,filterNum,ax,**kwargs):
        if not isinstance(filterNum,list):
            filterNum = [filterNum]
        filters = self.get(filterNum)
        for k in filters.keys():
            data = filters[k][1]
            if len(filterNum) == 1:
                ax.plot(data.wavelength,data.transmission,\
                            label=self.getID(int(k),short=False),**kwargs)
            else:
                ax.plot(data.wavelength,data.transmission,\
                            label=self.getID(int(k),short=False))
        return
            
