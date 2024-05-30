#!/usr/bin/python
# -*- coding: latin-1 -*-

## This two lines is to chose the econding
# =============================================================================
# Standard Python modules
# =============================================================================
import numpy as np
import h5py
import warnings
#import torch

###############################################################################

j = 1j

###############################################################################
## Parallelization that requieres mpi4py to be installed, if mpi4py was not installed successfully comment frome line 32 to line 40 (included)
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    size = comm.Get_size()
except:
    print('mpi4py is required for parallelization')
    myrank=0
    

# =============================================================================
# classes we will use
# =============================================================================

class globalparameters: 
    def __init__(self, t = None, freq = None, w = None):
        self.t = t
        self.freq = freq
        self.w = w

# =============================================================================

class inputdatafromfile: 
    def __init__(self, path, trace_start, trace_end, time_start, time_end, sample = 1):
        # sample = 1 if it is a sample and 0 if not and we want to compute the covariance with
        # trace_start, trace_end for the number of traces to study
        
        with h5py.File(path, "r") as f:
            self.time = np.array(f["timeaxis"])
            if trace_end == -1:
                self.numberOfTrace = len(f)-1 # on enleve timeaxis
                trace_end = len(f)-2
            else:
                self.numberOfTrace = trace_end-trace_start+1

            self.Pulseinit = [np.array(f[str(trace)]) for trace in range(trace_start, trace_end+1)] #list of all time trace without the ref, traces are ordered
            if sample:
                self.ref_number = self.choose_ref_number()
            else:
                self.ref_number = None
            self.timestamp = None
            try:
                self.timestamp = [f[str(trace)].attrs["TIMESTAMP"] for trace in range(trace_start, trace_end+1) ]
            except:
                pass
            
    def choose_ref_number(self):
        norm = []
        pseudo_norm = []
        temp = mean(self.Pulseinit)

        for i in range(self.numberOfTrace):
            pseudo_norm.append(np.dot(self.Pulseinit[i], temp))
            
        for i in range(self.numberOfTrace):
            norm.append(np.dot(self.Pulseinit[i], self.Pulseinit[i]))
        proximity = np.array(pseudo_norm)/np.array(norm)

        return np.abs(proximity - 1).argmin()  #proximite la plus proche de 1
    

class getreferencetrace:
    def __init__(self, path, ref_number, trace_start, time_start):
        with h5py.File(path, "r") as f:
            self.Pulseinit = np.array(f[str(trace_start+ref_number)])
            self.Spulseinit = torch_rfft(self.Pulseinit)  ## We compute the spectrum of the measured pulse
            
    def transferfunction(self, myinputdata):
        return self.Spulseinit/myinputdata.Spulse
        


# =============================================================================

class mydata:
     def __init__(self, pulse):
        self.pulse = pulse # pulse
        self.Spulse = torch_rfft((pulse)) # spectral field
        
        
class datalist:
     def __init__(self):
        self.pulse = []
        self.moyenne = []
        self.time_std = []
        self.freq_std = []
        self.freq_std_with_window = []
        self.freq_std_to_save = []  #if superresolution save the std of the cutted traces as output

        self.covariance = None
        self.covariance_inverse = None

        
     def add_trace(self, pulse):
        self.pulse.append(pulse) # pulse with sample
        
     def add_ref(self, ref_number, refpulse):
        self.pulse.insert(ref_number, refpulse) # pulse with sample
        

def torch_rfft(signal, axis = -1):
    #return torch.fft.rfft(torch.from_numpy(np.array(signal)), dim = axis).numpy()
    return np.fft.rfft(np.array(signal), axis = axis)
 

def torch_irfft(signal, axis = -1, n = None):
    #return torch.fft.irfft(torch.from_numpy(np.array(signal)), dim = axis, n = n).numpy()
    return np.fft.irfft(np.array(signal), axis = axis, n = n)


def mean(array):
    moyenne = np.zeros(len(array[-1]))
    for i in array:
        moyenne = moyenne+ i
    moyenne = moyenne/len(array)
    return moyenne

def std(array,moyenne):
    somme = np.zeros(len(array[-1]))
    for i in array:
        somme+= abs(i - moyenne)**2
    ecart = np.sqrt((somme/len(array)))
    return ecart

