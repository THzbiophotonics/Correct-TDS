#!/usr/bin/python
# -*- coding: latin-1 -*-

## This two lines is to chose the econding
# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
from Correct_TDS import fitf as TDS
import pickle
import subprocess
import numpy as np
import h5py
import warnings
import multiprocessing # For the break button
from scipy import signal
from sklearn.covariance import GraphicalLassoCV, LedoitWolf, OAS
from pathlib import Path as path_

import numba #pour inverser rapidement
@numba.jit
def inv_nla_jit(A):
    return np.linalg.inv(A)
###############################################################################

j = 1j

# =============================================================================
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

class ControlerBase:
    def __init__(self):
        self.clients_tab1 = list()
        self.clients_tab3 = list()
        self.message = ""


    def addClient(self, client):
        self.clients_tab1.append(client)

    def addClient3(self, client):
        self.clients_tab3.append(client)
        

    def refreshAll(self, message):
        self.message = message
        for client in self.clients_tab1:
            client.refresh()
            
    def refreshAll3(self, message):
        self.message = message
        for client in self.clients_tab3:
            client.refresh()


class Controler(ControlerBase):
    def __init__(self):
        super().__init__()
        # Initialisation:
        
        self.myreferencedata=TDS.getreferencetrace
        
        self.data=TDS.inputdatafromfile
        self.data_without_sample=TDS.inputdatafromfile
        
        self.myglobalparameters=TDS.globalparameters()
        
        self.nsample=None
        self.nsamplenotreal=None
        self.dt=None   ## Sampling rate
                
        self.myinput = TDS.datalist()
        self.mydatacorrection=TDS.datalist() 
        self.myinput_without_sample = TDS.datalist()
        
        self.ncm = None  #noise convolution matrix
        self.ncm_inverse = None
        
        self.delay_correction = []
        self.leftover_correction = []
        self.periodic_correction = []
        self.dilatation_correction = []
        
        ## parameters for the optimization algorithm
        self.swarmsize=1
        self.maxiter=1

        # Variables for existence of temp Files
        # self.is_temp_file_3 = 0 # temp file storing optimization results
        # self.is_temp_file_4 = 0
        # self.is_temp_file_5 = 0 # temp file storing algorithm choices
        

        # Variable to see if initialisation is done
        self.initialised = 0
        self.initialised_param = 0
        self.optim_succeed = 0
        
        self.path_data = None
        self.path_data_ref = None
        self.reference_number = None
        self.fit_delay = 0
        self.delaymax_guess = None
        self.delay_limit = None
        self.fit_leftover_noise = 0
        self.leftcoef_guess = None
        self.leftcoef_limit = None
        self.fit_periodic_sampling = 0
        self.periodic_sampling_freq_limit = None
        self.fit_dilatation = 0
        self.dilatation_limit = None
        self.dilatationmax_guess = None
        
        self.Lfiltering = None
        self.Hfiltering = None
        
        self.algo = None
        self.mode = None
        
        self.Freqwindow = None
        self.timeWindow = None
        self.fopt = []
        self.fopt_init = []
        
        self.trace_start = 0
        self.trace_end = -1
        self.time_start = 0
        self.time_end = -1
        
        self.optim=TDS.Optimization()
        self.myfitdata = TDS.myfitdata
        
        self.optimization_process=None
        self.interface_process=None
        
# =============================================================================
# Initialisation tab
# =============================================================================

    def init(self):
        self.refreshAll("Initialisation: Ok")

    def loading_text(self):
        self.refreshAll("\n Processing... \n")

    def choices_ini(self, path_data, path_data_ref, trace_start, trace_end, time_start, time_end,
                    Lfiltering_index, Hfiltering_index, cutstart, cutend,sharpcut, modesuper, apply_window):
        """Process all the informations given in the first panel of initialisation: 
            create instances of classes to store data, apply filters"""

        self.path_data = path_data
        self.path_data_ref = path_data_ref
        self.optim_succeed = 0
        
        self.myinput = TDS.datalist() # Warning: Don't remove the parentheses
        self.myinput_without_sample = TDS.datalist()
        self.mydatacorrection=TDS.datalist()
        
        self.ncm = None  
        self.ncm_inverse = None
        self.delay_correction = []
        self.leftover_correction = []
        self.periodic_correction = []
        self.dilatation_correction = []
        self.fopt = []
        self.fopt_init = []
        
        self.trace_start = trace_start
        self.trace_end = trace_end
        self.time_start = time_start
        self.time_end = time_end
        
         #####################################################################
         
        with h5py.File(path_data, "r") as f:
            l = len(f)-1
            if trace_start >= l or (trace_end!=-1 and trace_end >= l):
                self.refreshAll3("Invalid length, the file contains "+str(l)+" traces")        
        
        self.data= TDS.inputdatafromfile(path_data, self.trace_start,self.trace_end, self.time_start, self.time_end)    ## We load the signal of the measured pulse

        if path_data_ref:
            with h5py.File(path_data_ref, "r") as f:
                l = len(f)-1
                if trace_start >= l or (trace_end!=-1 and trace_end >= l):
                    self.refreshAll3("Invalid length, the file without sample contains "+str(l)+" traces")  
                    
            self.data_without_sample= TDS.inputdatafromfile(path_data_ref, self.trace_start,self.trace_end, self.time_start, self.time_end, sample = 0)
        
            if self.data_without_sample.numberOfTrace != self.data.numberOfTrace or len(self.data_without_sample.Pulseinit[0])!=len(self.data.Pulseinit[0]):
                self.refreshAll3("Invalid length, the two h5 files doesn't have the same dimension")
                
     ###########################################################################
    
    
        self.reference_number = self.data.ref_number
        self.myreferencedata = TDS.getreferencetrace(path_data, self.reference_number, self.trace_start, self.time_start)

        self.myglobalparameters.t = self.data.time*1e-12 # this assumes input files are in ps ## We load the list with the time of the experiment
        self.nsample = len(self.myglobalparameters.t)
        self.dt=self.myglobalparameters.t.item(2)-self.myglobalparameters.t.item(1)  ## Sample rate
        self.myglobalparameters.freq = np.fft.rfftfreq(self.nsample, self.dt)        ## We create a list with the frequencies for the spectrum
        self.myglobalparameters.w = self.myglobalparameters.freq*2*np.pi
        
        self.Lfiltering = Lfiltering_index 
        self.Hfiltering = Hfiltering_index 
        
        if modesuper == 1:
            frep=99.991499600e6 # repetition frequency of the pulse laser used in the tds measurments in Hz, 99
            nsampleZP=np.round(1/(frep*self.dt)) #number of time sample betwen two pulses. IT has to be noted that it could be better to have an integer number there then the rounding does not change much
            self.nsamplenotreal=nsampleZP.astype(int)
            self.myglobalparameters.t=np.arange(nsampleZP)*self.dt  # 0001 #
            self.myglobalparameters.freq = np.fft.rfftfreq(self.nsamplenotreal, self.dt)
            self.myglobalparameters.w = 2*np.pi*self.myglobalparameters.freq

        else:
            self.nsamplenotreal = self.nsample 
            
        Freqwindowstart = np.ones(len(self.myglobalparameters.freq))
        Freqwindowend = np.ones(len(self.myglobalparameters.freq))
        if self.Lfiltering:
            stepsmooth = cutstart/sharpcut
            Freqwindowstart = 0.5+0.5*np.tanh((self.myglobalparameters.freq-cutstart)/stepsmooth)
        if self.Hfiltering:
            stepsmooth = cutend/sharpcut
            Freqwindowend = 0.5-0.5*np.tanh((self.myglobalparameters.freq-cutend)/stepsmooth)                                 
        
        self.Freqwindow = Freqwindowstart*Freqwindowend
        self.timeWindow = np.ones(self.nsamplenotreal)
        
        for trace in  range(self.data.numberOfTrace):

            myinputdata=TDS.mydata(self.data.Pulseinit[trace])    ## We create a variable containing the data related to the measured pulse with sample

            if modesuper == 1:
                self.mode = "superresolution"
                myinputdata=TDS.mydata(np.pad(myinputdata.pulse,(0,self.nsamplenotreal-self.nsample),'constant',constant_values=(0)))

                if trace == 0: # on fait le padding une seule fois sur la ref
                    self.myreferencedata.Pulseinit=np.pad(self.myreferencedata.Pulseinit,(0,self.nsamplenotreal-self.nsample),'constant',constant_values=(0))
                    self.myreferencedata.Spulseinit=(TDS.torch_rfft((self.myreferencedata.Pulseinit)))    # fft computed with GPU
            else:
                self.mode = "basic"
                
            # Filter data
            if trace == 0:
                self.myreferencedata.Spulseinit = self.myreferencedata.Spulseinit*self.Freqwindow
                self.myreferencedata.Pulseinit  = TDS.torch_irfft(self.myreferencedata.Spulseinit, n = self.nsamplenotreal)

            myinputdata.Spulse         = myinputdata.Spulse        *self.Freqwindow
            myinputdata.pulse          = TDS.torch_irfft(myinputdata.Spulse, n = self.nsamplenotreal)
                    
                

            self.myinput.add_trace(myinputdata.pulse)
            if path_data_ref:
                self.myinput_without_sample.add_trace(self.data_without_sample.Pulseinit[trace])
                self.data_without_sample.Pulseinit[trace] = []
            
        if path_data_ref:
            self.myinput_without_sample.moyenne = np.mean(self.myinput_without_sample.pulse, axis= 0)
            
        self.myinput.moyenne = np.mean(self.myinput.pulse, axis= 0)  ### TDS.mean and TDS.std instead for big dataset
        self.myinput.time_std = np.std(self.myinput.pulse, axis = 0)
        self.myinput.freq_std = np.std(TDS.torch_rfft(self.myinput.pulse, axis = 1), axis = 0)
        
        if apply_window == 1:  # it's not a linear operation in freq domain
            windows = signal.tukey(self.nsamplenotreal, alpha = 0.05)
            self.myinput.freq_std_with_window = np.std(TDS.torch_rfft(self.myinput.pulse*windows, axis = 1), axis = 0)

        # if not os.path.isdir("temp"):
        #     os.mkdir("temp")
        # with open(os.path.join("temp",'temp_file_6.bin'),'wb') as f:
        #     pickle.dump(self.myinput,f,pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(self.myreferencedata,f,pickle.HIGHEST_PROTOCOL)
        
        # with open(os.path.join("temp",'temp_file_7.bin'),'wb') as f:
        #     pickle.dump(self.myglobalparameters,f,pickle.HIGHEST_PROTOCOL)
        #     pickle.dump(apply_window,f,pickle.HIGHEST_PROTOCOL)
        
        self.optim.vars_temp_file_6_data=self.myinput
        self.optim.vars_temp_file_6_ref=self.myreferencedata
        self.optim.vars_temp_file_7_globalparameters=self.myglobalparameters
        self.optim.vars_temp_file_7_apply_window=apply_window
        
        self.myfitdata.myglobalparameters = self.myglobalparameters
        
        self.data.Pulseinit = [] #don't forget to empty it, important for memory


    def choices_ini_param( self, fit_delay, delaymax_guess, delay_limit, fit_dilatation, dilatation_limit,dilatationmax_guess,
                    fit_periodic_sampling, periodic_sampling_freq_limit, fit_leftover_noise, leftcoef_guess, leftcoef_limit):


        self.fit_delay = fit_delay
        self.delaymax_guess = delaymax_guess
        self.delay_limit = delay_limit
        self.fit_leftover_noise = fit_leftover_noise
        self.leftcoef_guess = leftcoef_guess
        self.leftcoef_limit = leftcoef_limit
        self.fit_periodic_sampling = fit_periodic_sampling
        self.periodic_sampling_freq_limit = periodic_sampling_freq_limit
        self.fit_dilatation = fit_dilatation
        self.dilatation_limit = dilatation_limit
        self.dilatationmax_guess = dilatationmax_guess

            # files for choices made
        mode_choicies_opt=[self.path_data, self.path_data_ref, self.reference_number, self.fit_dilatation, self.dilatation_limit, self.dilatationmax_guess,
                               self.Freqwindow,self.timeWindow, self.fit_delay, self.delaymax_guess, self.delay_limit,  self.mode, self.nsample,
                               self.fit_periodic_sampling, self.periodic_sampling_freq_limit, self.fit_leftover_noise, self.leftcoef_guess, self.leftcoef_limit]
    
        # # if not os.path.isdir("temp"):
        # #     os.mkdir("temp")
        # with open(os.path.join("temp",'temp_file_1_ini.bin'),'wb') as f:
        #     pickle.dump(mode_choicies_opt,f,pickle.HIGHEST_PROTOCOL)
        self.optim.vars_temp_file_1_ini=mode_choicies_opt
        self.myfitdata.vars_temp_file_1_ini=mode_choicies_opt




# =============================================================================
# Model parameters
# =============================================================================
   
    def invalid_swarmsize(self):
        self.refreshAll3("Invalid swarmsize. \n")

    def invalid_niter(self):
        self.refreshAll3("Invalid number of iterations. \n")
        
        
# =============================================================================
# Optimization
# =============================================================================
        
    def algo_parameters(self,choix_algo,swarmsize,niter, niter_ps):
        """Save algorithm choices in temp file 5"""
        self.algo=choix_algo
        mode_choicies_opt=[choix_algo,int(swarmsize),int(niter), int(niter_ps)]
        # if not os.path.isdir("temp"):
        #     os.mkdir("temp")

        # with open(os.path.join("temp",'temp_file_5.bin'),'wb') as f:
        #     pickle.dump(mode_choicies_opt,f,pickle.HIGHEST_PROTOCOL)
        # self.is_temp_file_5 = 1
        self.optim.vars_temp_file_5=mode_choicies_opt

        self.refreshAll3("")
        
        
    def begin_optimization(self,nb_proc):
        
        
        self.ncm = None
        
        """Run optimization and update layers"""
        output=""
        error=""
        returncode=0
        self.optim.interrupt=False
        self.optim.optimize(nb_proc)  
   
   
        # #Creating an optimisation process
        # self.optimization_process = multiprocessing.Process(target=self.optim.optimize, args=(nb_proc,),daemon=True)
        # self.optimization_process.start()
        
        # #waiting for the end of the optimisation process
        # self.optimization_process.join()
        # print("Optimisation terminée.")

        
#         if sys.platform=="win32" or sys.platform=="cygwin":
#             print("OS:Windows \n")
#             if not os.path.isdir("temp"):
#                 os.mkdir("temp")
#             optimization_filename = os.path.join('temp',"opt.bat")
#             try:
#                 with open(optimization_filename, 'w') as OPATH:
#                     OPATH.writelines(['call set Path=%Path%;C:\ProgramData\Anaconda3 \n',
#                     'call set Path=%Path%;C:\ProgramData\Anaconda3\condabin \n',
#                     'call set Path=%Path%;C:\ProgramData\Anaconda3\Scripts \n',
#                     #'call conda activate \n', 
#                     f'call mpiexec -n {nb_proc} python opt.py'])
# #                    OPATH.writelines([f'call mpiexec -n {nb_proc} opt.exe'])
#                 subprocess.call(optimization_filename)
#                 returncode = 0
#                 error = ""
#                 output = ""
#             except:
#                 print("No parallelization! You don't have MPI installed or there's a problem with your MPI.")
#                 with open(optimization_filename, 'w') as OPATH:
#                     OPATH.writelines([f'call opt.exe'])
#                 subprocess.call(optimization_filename)
#                 returncode = 0
#                 error = ""
#                 output = ""
#         elif sys.platform=="linux" or sys.platform=="darwin":
#             # print("OS:Linux/MacOS \n")
#             optimization_filename = os.path.join('temp',"opt.sh")
#             # print(optimization_filename)
#             try:
#                 # Check if Open MPI is correctly installed
#                 try:
#                     command = 'mpiexec --version'
#                     process=subprocess.Popen(command.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#                     output_mpi,error_mpi = process.communicate()
#                     returncode_mpi=process.returncode
#                 except:
#                     returncode_mpi = 1
#                     error_mpi = "Command mpiexec not recognized."

#                 try:
#                     command = './correct-env/bin/python --version'
#                     process=subprocess.Popen(command.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#                     output_py3,error_py3 = process.communicate()
#                     returncode_py3=process.returncode
#                     # print(f"####\npython venv\n####")
#                     python_path = "./correct-env/bin/python"
#                 except:
#                     try:
#                         command = "python --version"
#                         process=subprocess.Popen(command.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#                         output_py3,error_py3 = process.communicate()
#                         returncode_py3=process.returncode
#                         python_path = "python"
#                         # print(f" python version worked" if returncode_py3 == 0 else f"python version failed")
#                     except:
#                         returncode_py3 = 1
#                         error_py3 = "Command python not recognized."

#                 # Run optimization
#                 if returncode_mpi==0:
#                     if returncode_py3==0:
#                         # fpath = os.getcwd()
#                         # fpath = os.path.join(fpath,"opt.py")
#                         # print(fpath)
#                         # print(f"enter mpiexec command")
#                         # command = f"#!/bin/sh\nmpiexec -n 1 {python_path} opt.py"
#                         command = f"#!/bin/sh\npython opt.py"
#                         # command = 'mpiexec -n 20 {1} opt.py'.format(python_path)
#                         # print(command)
#                     else:
#                         print(f"Problem with python command : \n {error_py3} \n")
#                         return(0)
#                 else:
#                     print(f"No parallelization! You don't have MPI installed or there's a problem with your MPI: \n {error_mpi}")
#                     if returncode_py3==0:
#                         command = f'#!/bin/sh\n{python_path} opt.py'
#                     else:
#                         print(f"Problem with python command : \n {error_py3} \n")
#                         return(0)

#                 try:
#                     with open(optimization_filename, 'w') as OPATH:
#                         OPATH.writelines(command)
#                     # print(f"permission ask")
#                     returncode = subprocess.call(f'chmod +x ./{optimization_filename}',shell=True)
#                     if returncode == 0:
#                         # print(f"permission granted")
#                         # print(os.getcwd())
#                         # returncode = subprocess.call('./{}'.format(optimization_filename),shell=True)
#                         optimization_filename = path_().cwd().joinpath("temp").joinpath("opt.sh")
#                         # print(f"{optimization_filename}")
#                         process = subprocess.run(f"{optimization_filename}",shell=True)
#                         # print(process.returncode)
#                         # print(returncode)
#                     if returncode == 1:
#                         command = ""
#                         with open("launch_opt.py", 'w') as OPATH:
#                             OPATH.writelines(command)
#                             # print(f"created launch_opt.py")
#                         try:
#                             import launch_opt
#                             try:
#                                 f=open(os.path.join("temp",'temp_file_3.bin'),'rb')
#                                 f.close()
#                                 returncode=0
#                             except:
#                                 print("Unknown problem. cannot open temp file 3")
#                                 sys.exit()
#                         except:
#                             print("Unknown problem. cannot import launch opt")
#                             sys.exit()
#                 except:
#                     returncode = 1
#                     error = "Unknow problem."
#                     output = ""
#             except:
#                 print("Unknown problem.")
#                 sys.exit()

#         else:
#             print("System not supported.")
#             return(0)

#         # We dont verify if returncode==0 anymore
#         if returncode==0:
#         with open(os.path.join("temp",'temp_file_3.bin'),'rb') as f:
        
        # Vérification de vars_temp_file_3
        # var_inter = self.optim.vars_temp_file_3
        # print(f"var_inter before looks like : {var_inter}")
        # f=open(os.path.join("temp",'temp_file_3.bin'),'rb')
        self.is_temp_file_3 = 1
        self.delay_correction = []
        self.leftover_correction = []
        self.periodic_correction = []
        self.dilatation_correction = []
        self.fopt = []
        sum_fopt = 0
        if self.fit_delay or self.fit_leftover_noise or self.fit_dilatation:
            for i in range(self.data.numberOfTrace): # we begin at 2 cause the first values don't work
                var_inter=self.optim.vars_temp_file_3
                temp_var_inter = var_inter[i]
                xopt=temp_var_inter[0] # replace var_inter w temp_var_inter
                fopt=temp_var_inter[1] # replace var_inter w temp_var_inter
                self.fopt.append(fopt)
                if self.fit_leftover_noise:
                    self.leftover_correction.append(xopt[-2:])
                if self.fit_delay:
                    self.delay_correction.append(xopt[0])

                if self.fit_dilatation:
                    if self.fit_delay:
                        self.dilatation_correction.append(xopt[1:3])
                    else:
                        self.dilatation_correction.append(xopt[0:2])
                        
            sum_fopt = np.sum(self.fopt)
        
        if self.fit_periodic_sampling:
            var_inter=self.optim.vars_temp_file_3
            xopt_ps=var_inter[0]
            fopt_ps=var_inter[1]
            self.periodic_correction.append(xopt_ps)
            
        
        # with open(os.path.join("temp",'temp_file_2.bin'),'rb') as f:
        self.mydatacorrection = self.optim.vars_temp_file_2_datacorrection
        self.fopt_init = self.optim.vars_temp_file_2_fopt       #available only for delay,dilatation,amplitude correction
        

        message = "Optimization terminated successfully\nCheck the output directory for the result\n\n"
  
        if self.fit_periodic_sampling:
        #     message += f'For periodic sampling: \n The best error was:     {fopt_ps}' + f'\nThe best parameters were:     {xopt_ps}\n' + "\n"
        # if self.fit_leftover_noise or self.fit_delay or self.fit_dilatation:
            message+= f"Sum of std Pulse E field (sqrt(sum(std^2))):\n   before correction \t{np.sqrt(sum(self.myinput.time_std**2))}\n"
            message+= f"   after correction \t{np.sqrt(sum(self.mydatacorrection.time_std**2))}\n\n"
            
            #message+= "Sum of errors :\n   before correction \t{}\n".format(sum_fopt)
            #message+= " Sum of errors after correction \t{}\n\n".format(sum_fopt)

        citation= "Please cite this paper in any communication about any use of Correct@TDS : \nComing soon..."
        message += citation
        
        self.refreshAll3(message)
        # else:
        #     self.refreshAll3(f"Output : \n {output} \n")
        #     print("System not supported. \n")
        #     print(f'Output : \n {output} \n Error : \n {error} \n')
        #     return(0)
        
    def stop_optimization(self):
        self.optim.interrupt=True
        # if self.optimization_process is not None:
        #     self.optimization_process.join()
        
    def save_data(self, filename, path, file, cov_algo = 3):
      
        citation= "Please cite this paper in any communication about any use of Correct@TDS : \n Coming soon..."
        custom = "\n Average over "+str(self.data.numberOfTrace)+" waveforms. Timestamp: "
        try:
            if self.initialised: 
                    if self.optim_succeed:
                        if file == 0:
                            title = "\n timeaxis (ps) \t E-field"
                            out = np.column_stack((self.data.time, self.mydatacorrection.moyenne[:self.nsample]))

                            if self.data.timestamp:
                                custom+= str(self.data.timestamp[0])
                            else:
                                custom+= "unknown"
                            np.savetxt(os.path.join(path,filename),out, header= citation+custom+title, delimiter = "\t")
                    
                        
                        if file == 1:
                            to_save = []
                            if self.fit_delay:
                                to_save.append(self.delay_correction)
                            else:
                                to_save.append([0]*self.data.numberOfTrace) #empty
                                
                            if self.fit_dilatation:
                                to_save.append([self.dilatation_correction[i][0] for i in range(self.data.numberOfTrace)])
                                to_save.append([self.dilatation_correction[i][1] for i in range(self.data.numberOfTrace)])
                            else:
                                to_save.append([0]*self.data.numberOfTrace) #empty
                                to_save.append([0]*self.data.numberOfTrace) #empty
    
                                
                            if self.fit_leftover_noise:
                                to_save.append([self.leftover_correction[i][0] for i in range(self.data.numberOfTrace)])
                            else:
                                to_save.append([0]*self.data.numberOfTrace) #empty
                                
                            if self.fit_periodic_sampling:
                                to_save.append([self.periodic_correction[0][0]]*self.data.numberOfTrace)  
                                to_save.append([self.periodic_correction[0][1]]*self.data.numberOfTrace)  
                                to_save.append([self.periodic_correction[0][2]]*self.data.numberOfTrace)  
                            else:
                                to_save.append([0]*self.data.numberOfTrace) #empty
                                to_save.append([0]*self.data.numberOfTrace) #empty
                                to_save.append([0]*self.data.numberOfTrace) #empty
    
                            np.savetxt(os.path.join(path,filename), np.transpose(to_save), delimiter = "\t", header=citation+"\n delay \t alpha_dilatation \t beta_dilatation \t a_amplitude \t A \t nu \t phi")
                
                        if file == 2:
                            #save in hdf5
                            hdf =  h5py.File(os.path.join(path,filename),"w")
                            dataset = hdf.create_dataset(str("0"), data = self.mydatacorrection.pulse[0][:self.nsample])
                            dataset.attrs["CITATION"] = citation
                            
                            if self.data.timestamp:
                                dataset.attrs["TIMESTAMP"] = self.data.timestamp[0]

                            hdf.create_dataset('timeaxis', data = self.data.time)
                                
                            count = 1
                            for i in self.mydatacorrection.pulse[1:]:
                                dataset = hdf.create_dataset(str(count), data = i[:self.nsample])

                                dataset.attrs["CITATION"] = citation
                                if self.data.timestamp:
                                    dataset.attrs["TIMESTAMP"] = self.data.timestamp[count]
                                count+=1
                            hdf.close()
                            
                        if file == 3:
                            title = "\n timeaxis (ps) \t Std E-field"

                            out = np.column_stack((self.data.time, self.mydatacorrection.time_std[:self.nsample]))
                            
                            if self.data.timestamp:
                                custom+= str(self.data.timestamp[0])
                            else:
                                custom+= "unknown"
                            np.savetxt(os.path.join(path,filename),out, header= citation+custom+title, delimiter = "\t")
                            
                        if file == 4:
                            title = "\n Frequency (Hz) \t Std E-field"
                            if self.mode == "superresolution":
                                if not self.mydatacorrection.freq_std_to_save:
                                    self.mydatacorrection.freq_std_to_save = np.std(TDS.torch_rfft([self.mydatacorrection.pulse[i][:self.nsample] for i in range(self.data.numberOfTrace)], axis = 1),axis = 0)
                                out = np.column_stack((np.fft.rfftfreq(self.nsample, self.dt), self.mydatacorrection.freq_std_to_save))
                            else:
                                out = np.column_stack((self.myglobalparameters.freq, self.mydatacorrection.freq_std))
                            
                            if self.data.timestamp:
                                custom+= str(self.data.timestamp[0])
                            else:
                                custom+= "unknown"
                            np.savetxt(os.path.join(path,filename),out, header= citation+custom+title, delimiter = "\t")
                            
                            
                        if file == 5:
                            if self.mydatacorrection.covariance is None:
                                if self.path_data_ref:
                                    transfer_function = TDS.torch_irfft(TDS.torch_rfft(self.mydatacorrection.moyenne[:self.nsample])/TDS.torch_rfft(self.myinput_without_sample.moyenne))
                            
                                if cov_algo == 1:
                                    if self.path_data_ref:
                                        cov_with_ref = LedoitWolf().fit([np.convolve(transfer_function, self.myinput_without_sample.pulse[i])[:self.nsample] for i in range(self.data.numberOfTrace) ])
                                        self.myinput_without_sample.covariance =  cov_with_ref.covariance_ /self.data.numberOfTrace
                                        cov_with_ref = []
                                        
                                    cov = LedoitWolf().fit(np.array(self.mydatacorrection.pulse)[:,:self.nsample])
                                    self.mydatacorrection.covariance =  cov.covariance_ /self.data.numberOfTrace
                                    

                                    
                                elif cov_algo == 2:
                                    if self.path_data_ref:
                                        cov_with_ref = OAS().fit([np.convolve(transfer_function, self.myinput_without_sample.pulse[i])[:self.nsample] for i in range(self.data.numberOfTrace) ])
                                        self.myinput_without_sample.covariance =  cov_with_ref.covariance_ /self.data.numberOfTrace
                                        cov_with_ref = []
                                        
                                    cov = OAS().fit(np.array(self.mydatacorrection.pulse)[:,:self.nsample])
                                    self.mydatacorrection.covariance = cov.covariance_ /self.data.numberOfTrace
                                    

                                    
                                elif cov_algo == 3:
                                    if self.path_data_ref:
                                        model = GraphicalLassoCV(cv = 3, alphas=2, n_refinements=10, max_iter = 100, mode = "cd", n_jobs=1, tol = 1e-4, verbose = True)
                                        cov_with_ref = model.fit([np.convolve(transfer_function, self.myinput_without_sample.pulse[i])[:self.nsample] for i in range(self.data.numberOfTrace) ])
                                        self.myinput_without_sample.covariance =  cov_with_ref.covariance_ /self.data.numberOfTrace
                                        alpha_myinput_without_sample = cov_with_ref.alpha_
                                        cov_with_ref = []
                                        
                                    model = GraphicalLassoCV(cv = 3, alphas=2, n_refinements=10, max_iter = 100, mode = "cd", n_jobs=1, tol = 1e-4, verbose = True)
                                    cov = model.fit(np.array(self.mydatacorrection.pulse)[:,:self.nsample])
                                    self.mydatacorrection.covariance = cov.covariance_ /self.data.numberOfTrace
                                    alpha_mydatacorrection = cov.alpha_
                                    
                                   
                            if self.path_data_ref:
                                hdf =  h5py.File(os.path.join(path,filename),"w")
                                self.ncm = self.mydatacorrection.covariance + self.myinput_without_sample.covariance
                                self.ncm_inverse = inv_nla_jit(self.ncm)
                                dataset = hdf.create_dataset("noise_convolution_matrix_inverse", data = self.ncm_inverse)
                                dataset.attrs["CITATION"] = citation
                                if self.data.timestamp:
                                    dataset.attrs["REFEFERENCE_TIMESTAMP"] = self.data.timestamp[0]
                                if self.data_without_sample.timestamp:
                                    dataset.attrs["SAMPLE_TIMESTAMP"] = self.data_without_sample.timestamp[0]
                                    
                            else:
                                hdf =  h5py.File(os.path.join(path,filename),"w")
                                if self.mydatacorrection.covariance_inverse is None:
                                    self.mydatacorrection.covariance_inverse = cov.precision_*self.data.numberOfTrace
                                    cov = []
                                dataset = hdf.create_dataset("covariance_inverse", data = self.mydatacorrection.covariance_inverse)
                                dataset.attrs["CITATION"] = citation
                                if self.data.timestamp:
                                    dataset.attrs["TIMESTAMP"] = self.data.timestamp[0]
                            hdf.close()

                    else: #if no optimization
                        if file == 0:
                            title = "\n timeaxis (ps) \t E-field"
                            if self.data.timestamp:
                                custom+= str(self.data.timestamp[0])
                            else:
                                custom+= "unknown"
                            out = np.column_stack((self.data.time, self.myinput.moyenne[:self.nsample]))

                            np.savetxt(os.path.join(path,filename),out, delimiter = "\t", header= citation+custom+title)                    
                        
                        if file == 2:
                            #save in hdf5
                            hdf =  h5py.File(os.path.join(path,filename),"w")
                            dataset = hdf.create_dataset(str("0"), data = self.myinput.pulse[0][:self.nsample])
                            dataset.attrs["CITATION"] = citation
                            if self.data.timestamp:
                                dataset.attrs["TIMESTAMP"] = self.data.timestamp[0]
                                
                            hdf.create_dataset('timeaxis', data = self.data.time)
                  
                            count = 1
                            for i in self.myinput.pulse[1:]:
                                dataset = hdf.create_dataset(str(count), data = i[:self.nsample])

                                dataset.attrs["CITATION"] = citation
                                if self.data.timestamp:
                                    dataset.attrs["TIMESTAMP"] = self.data.timestamp[count]
                                count+=1
                            hdf.close()
                            
                        if file == 3:
                            title = "\n timeaxis (ps) \t Std E-field"
                            if self.data.timestamp:
                                custom+= str(self.data.timestamp[0])
                            else:
                                custom+= "unknown"

                            out = np.column_stack((self.data.time, self.myinput.time_std[:self.nsample]))

                            np.savetxt(os.path.join(path,filename),out, delimiter = "\t", header= citation+custom+title)
                            
                            
                        if file == 4:
                            title = "\n Frequency (Hz) \t Std E-field"
                            if self.data.timestamp:
                                custom+= str(self.data.timestamp[0])
                            else:
                                custom+= "unknown"
                            if self.mode == "superresolution":
                                if not self.myinput.freq_std_to_save:
                                    self.myinput.freq_std_to_save = np.std(TDS.torch_rfft([self.myinput.pulse[i][:self.nsample] for i in range(self.data.numberOfTrace)], axis = 1),axis = 0)
                                out = np.column_stack((np.fft.rfftfreq(self.nsample, self.dt), self.myinput.freq_std_to_save))
                            else:
                                out = np.column_stack((self.myglobalparameters.freq, self.myinput.freq_std))

                            np.savetxt(os.path.join(path,filename),out, delimiter = "\t", header= citation+custom+title)
                            
                        if file == 5:
                            if self.myinput.covariance is None:
                                if self.path_data_ref:
                                    transfer_function = TDS.torch_irfft(TDS.torch_rfft(self.myinput.moyenne[:self.nsample])/TDS.torch_rfft(self.myinput_without_sample.moyenne))
                            
                                if cov_algo == 1:
                                    if self.path_data_ref:
                                        cov_with_ref = LedoitWolf().fit([np.convolve(transfer_function, self.myinput_without_sample.pulse[i])[:self.nsample] for i in range(self.data.numberOfTrace) ])
                                        self.myinput_without_sample.covariance =  cov_with_ref.covariance_ /self.data.numberOfTrace
                                        cov_with_ref = []
                                        
                                    cov = LedoitWolf().fit(np.array(self.myinput.pulse)[:,:self.nsample])
                                    self.myinput.covariance =  cov.covariance_ /self.data.numberOfTrace
                                    

                                    
                                elif cov_algo == 2:
                                    if self.path_data_ref:
                                        cov_with_ref = OAS().fit([np.convolve(transfer_function, self.myinput_without_sample.pulse[i])[:self.nsample] for i in range(self.data.numberOfTrace) ])
                                        self.myinput_without_sample.covariance =  cov_with_ref.covariance_ /self.data.numberOfTrace
                                        cov_with_ref = []
                                        
                                    cov = OAS().fit(np.array(self.myinput.pulse)[:,:self.nsample])
                                    self.myinput.covariance = cov.covariance_ /self.data.numberOfTrace
                                    

                                    
                                elif cov_algo == 3:
                                    if self.path_data_ref:
                                        model = GraphicalLassoCV(cv = 3, alphas=2, n_refinements=10, max_iter = 100, mode = "cd", n_jobs=1, tol = 1e-4, verbose = True)
                                        cov_with_ref = model.fit([np.convolve(transfer_function, self.myinput_without_sample.pulse[i])[:self.nsample] for i in range(self.data.numberOfTrace) ])
                                        self.myinput_without_sample.covariance =  cov_with_ref.covariance_ /self.data.numberOfTrace
                                        alpha_myinput_without_sample = cov_with_ref.alpha_
                                        cov_with_ref = []
                                        
                                    model = GraphicalLassoCV(cv = 3, alphas=2, n_refinements=10, max_iter = 100, mode = "cd", n_jobs=1, tol = 1e-4, verbose = True)
                                    cov = model.fit(np.array(self.myinput.pulse)[:,:self.nsample])
                                    self.myinput.covariance = cov.covariance_ /self.data.numberOfTrace
                                    alpha_myinput = cov.alpha_
                                    
                                   
                            if self.path_data_ref:
                                hdf =  h5py.File(os.path.join(path,filename),"w")
                                self.ncm = self.myinput.covariance + self.myinput_without_sample.covariance
                                self.ncm_inverse = inv_nla_jit(self.ncm)
                                dataset = hdf.create_dataset("noise_convolution_matrix_inverse", data = self.ncm_inverse)
                                dataset.attrs["CITATION"] = citation
                                if self.data.timestamp:
                                    dataset.attrs["REFEFERENCE_TIMESTAMP"] = self.data.timestamp[0]
                                if self.data_without_sample.timestamp:
                                    dataset.attrs["SAMPLE_TIMESTAMP"] = self.data_without_sample.timestamp[0]
                                    
                            else:
                                hdf =  h5py.File(os.path.join(path,filename),"w")
                                if self.myinput.covariance_inverse is None:
                                    self.myinput.covariance_inverse = cov.precision_*self.data.numberOfTrace
                                    cov = []
                                dataset = hdf.create_dataset("covariance_inverse", data = self.myinput.covariance_inverse)
                                dataset.attrs["CITATION"] = citation
                                if self.data.timestamp:
                                    dataset.attrs["TIMESTAMP"] = self.data.timestamp[0]
                            hdf.close()
                    return 1
            else:
                self.refreshAll3("Please enter initialization data first")
                return 0
        except Exception as e:
            self.refreshAll3("Please enter initialization data first")
            print(e)
            return 0

            
    def loading_text3(self):
        self.refreshAll3("\n Processing... \n")

    def message_log_tab3(self,message):
        self.refreshAll3(message)

    def error_message_path3(self):
        self.refreshAll3("Please enter a valid file(s). This file is maybe empty")

    def error_message_output_paths(self):
        self.refreshAll3("Invalid output paths.")

    def error_message_output_filename(self):
        self.refreshAll3("Invalid output filename.")


    def ploting_text3(self,message):
        self.refreshAll3(message)
    
    
    def no_temp_file_5(self):
        self.refreshAll3("Unable to execute without optimization parameters")
        
        
    def loading_text_break(self):
        self.refreshAll3("\n Breaking... \n")



