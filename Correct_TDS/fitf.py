#!/usr/bin/python
# -*- coding: latin-1 -*-

## This two lines is to chose the econding
# =============================================================================
# Standard Python modules
# =============================================================================
import numpy as np
import h5py
import warnings
import os, time
import pickle
from pyswarm import pso   ## Library for optimization
import scipy.optimize as optimize  ## Library for optimization
from scipy import signal
import json

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
# useful functions
# =============================================================================

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
        
        
class myfitdata: 
    def __init__(self, myinputdata, x):
        self.pulse = self.fit_input(myinputdata, x)
        self.Spulse = (torch_rfft(self.pulse))
        
        self.myreferencedata=None
        self.minval=None
        self.maxval=None
        self.myinputdata=None
        self.myglobalparameters=None
        self.dt=None
        
        self.vars_temp_file_1_ini = None
        
        
    def fit_input(self, myinputdata, x):
        # global mode, myglobalparameters, fit_delay, fit_leftover_noise, fit_dilatation, dt
        # with open(os.path.join("temp",'temp_file_1_ini.bin'),'rb') as f:
        #     [path_data, path_data_ref, reference_number, fit_dilatation, dilatation_limit, dilatationmax_guess, 
        #       freqWindow, timeWindow, fit_delay, delaymax_guess, delay_limit, mode, nsample,
        #       fit_periodic_sampling, periodic_sampling_freq_limit, fit_leftover_noise, leftcoef_guess, leftcoef_limit]=pickle.load(f)
            
        # with open(os.path.join("temp",'temp_file_7.bin'),'rb') as f:
        #     self.myglobalparameters = pickle.load(f)
        # self.myglobalparameters = self.vars_temp_file_7_globalparameters
        
        fit_dilatation=self.vars_temp_file_1_ini[3]
        fit_leftover_noise=self.vars_temp_file_1_ini[15]
        fit_delay=self.vars_temp_file_1_ini[8]
        
        dt=self.myglobalparameters.t.item(2)-self.myglobalparameters.t.item(1)   ## Sample rate


        
        leftover_guess = np.zeros(2)
        delay_guess = 0
        dilatation_coefguess = np.zeros(2)
        
        coef = np.zeros(2) #[a,c]
                
        if fit_delay:
            delay_guess = x[0]
            
        if fit_dilatation:
            if fit_delay:
                dilatation_coefguess = x[1:3]
            else:
                dilatation_coefguess = x[0:2]
            dilatation_coefguess[1] = 0
        if fit_leftover_noise:
                leftover_guess = x[-2:]

        coef[0] = leftover_guess[0] #a
        coef[1] = leftover_guess[1] #c

        Z = np.exp(j*self.myglobalparameters.w*delay_guess)
        myinputdatacorrected_withdelay = torch_irfft(Z*myinputdata.Spulse, n = len(self.myglobalparameters.t))

        leftnoise = np.ones(len(self.myglobalparameters.t)) - coef[0]*np.ones(len(self.myglobalparameters.t))   #(1-a)    
        myinputdatacorrected = leftnoise*(myinputdatacorrected_withdelay  
                                          - (dilatation_coefguess[0]*self.myglobalparameters.t)*np.gradient(myinputdatacorrected_withdelay, dt))
        

        return myinputdatacorrected




        
        
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
        






#### Classes for callbacks

class Callback_bfgs(object):
    def __init__(self,monerreur):
        self.nit = 0
        self.monerreur=monerreur
        
    def __call__(self, par, convergence=0):
        self.nit += 1
        with open('algo_bfgs_out.txt', 'a+') as filehandle:
            filehandle.write('\n iteration number %d ; error %s ; parameters %s \r\n' % (self.nit, self.monerreur(par), par))
            
class Callback_slsqp(object):
    def __init__(self,monerreur):
        self.nit = 0
        self.monerreur=monerreur
        
    def __call__(self, par, convergence=0):
        self.nit += 1
        with open('algo_slsqp_out.txt', 'a+') as filehandle:
            filehandle.write('\n iteration number %d ; error %s ; parameters %s \r\n' % (self.nit, self.monerreur(par), par))


class Callback_annealing(object):
    def __init__(self,monerreur):
        self.nit = 0
        self.monerreur=monerreur
        
    def __call__(self, par, f, context):
        self.nit += 1
        with open('algo_dualannealing_out.txt', 'a+') as filehandle:
            filehandle.write('\n iteration number %d ; error %s ; parameters %s \r\n' % (self.nit, self.monerreur(par), par))




class Optimization():
    def __init__(self):
        
        #Stocking all the variables which were in temp files, in majority used to optimize
        self.vars_temp_file_1_ini = None
        self.vars_temp_file_2_datacorrection = None
        self.vars_temp_file_2_fopt = None
        self.vars_temp_file_3 = None
        self.vars_temp_file_5 = None
        self.vars_temp_file_6_data = None
        self.vars_temp_file_6_ref = None
        self.vars_temp_file_7_globalparameters = None
        self.vars_temp_file_7_apply_window = None
        
        self.myreferencedata=None
        self.minval=None
        self.maxval=None
        self.myinputdata=None
        self.myglobalparameters=None
        self.dt=None
        self.mymean=None
        
        self.interrupt=False #Becoming True when break button pressed
        
        
    def errorchoice(self):
        # global mode, myglobalparameters, myinputdata, fit_delay, fit_leftover_noise, fit_dilatation, dt, nsample, maxval, minval
        nsample=self.vars_temp_file_1_ini[12]
        fit_dilatation=self.vars_temp_file_1_ini[3]
        fit_leftover_noise=self.vars_temp_file_1_ini[15]
        fit_delay=self.vars_temp_file_1_ini[8]
        
        def monerreur(x):
            
            leftover_guess = np.zeros(2)
            delay_guess = 0
            dilatation_coefguess = np.zeros(2)
            
            coef = np.zeros(2) #[a,c]

            x = x*(self.maxval-self.minval)+self.minval
            
            if fit_delay:
                delay_guess = x[0]
                
            if fit_dilatation:
                if fit_delay:
                    dilatation_coefguess = x[1:3]
                else:
                    dilatation_coefguess = x[0:2]
                dilatation_coefguess[1] = 0

            if fit_leftover_noise:
                    leftover_guess = x[-2:]

            coef[0] = leftover_guess[0] #a
            coef[1] = leftover_guess[1] #c

            Z = np.exp(j*self.myglobalparameters.w*delay_guess)
            myinputdatacorrected_withdelay = torch_irfft(Z*self.myinputdata.Spulse, n = len(self.myglobalparameters.t))

            leftnoise = np.ones(len(self.myglobalparameters.t)) - coef[0]*np.ones(len(self.myglobalparameters.t))   #(1-a)    
            myinputdatacorrected = leftnoise*(myinputdatacorrected_withdelay 
                                              - (dilatation_coefguess[0]*self.myglobalparameters.t)*np.gradient(myinputdatacorrected_withdelay, self.dt))
            erreur = np.linalg.norm(self.myreferencedata.Pulseinit[:nsample] - myinputdatacorrected[:nsample] )/np.linalg.norm(self.myreferencedata.Pulseinit[:nsample])
            return erreur
                
        return monerreur
        


    def fit_input(self, myinputdata, x):
        # global mode, myglobalparameters, fit_delay, fit_leftover_noise, fit_dilatation, dt
        
        fit_dilatation=self.vars_temp_file_1_ini[3]
        fit_leftover_noise=self.vars_temp_file_1_ini[15]
        fit_delay=self.vars_temp_file_1_ini[8]
        
        leftover_guess = np.zeros(2)
        delay_guess = 0
        dilatation_coefguess = np.zeros(2)
        
        coef = np.zeros(2) #[a,c]
                
        if fit_delay:
            delay_guess = x[0]
            
        if fit_dilatation:
            if fit_delay:
                dilatation_coefguess = x[1:3]
            else:
                dilatation_coefguess = x[0:2]
            dilatation_coefguess[1] = 0
        if fit_leftover_noise:
                leftover_guess = x[-2:]

        coef[0] = leftover_guess[0] #a
        coef[1] = leftover_guess[1] #c

        Z = np.exp(j*self.myglobalparameters.w*delay_guess)
        myinputdatacorrected_withdelay = torch_irfft(Z*self.myinputdata.Spulse, n = len(self.myglobalparameters.t))

        leftnoise = np.ones(len(self.myglobalparameters.t)) - coef[0]*np.ones(len(self.myglobalparameters.t))   #(1-a)    
        myinputdatacorrected = leftnoise*(myinputdatacorrected_withdelay  
                                          - (dilatation_coefguess[0]*self.myglobalparameters.t)*np.gradient(myinputdatacorrected_withdelay, self.dt))

        return myinputdatacorrected
    
    # =============================================================================
    # def errorchoice_pyOpt(): 
    #     def objfunc(x):  ## Function used in the Optimization function from pyOpt. For more details see http://www.pyopt.org/quickguide/quickguide.html
    #         optim = Optimization()
    #         monerreur = optim.errorchoice
    #         f = monerreur(x)
    #         fail = 0
    #         return f, 1, fail
    #     return objfunc


    # =============================================================================
        
    def optimize(self,nb_proc):
        # =============================================================================
        # We load the model choices
        # =============================================================================
        [path_data, path_data_ref, reference_number, fit_dilatation, dilatation_limit, dilatationmax_guess, 
         freqWindow, timeWindow, fit_delay, delaymax_guess, delay_limit, mode, nsample,
         fit_periodic_sampling, periodic_sampling_freq_limit, fit_leftover_noise, leftcoef_guess, leftcoef_limit]=self.vars_temp_file_1_ini

        # data = datalist()
        
        data = self.vars_temp_file_6_data
        
        self.myreferencedata=self.vars_temp_file_6_ref

        self.myglobalparameters = globalparameters()
        
        self.myglobalparameters = self.vars_temp_file_7_globalparameters
        apply_window = self.vars_temp_file_7_apply_window
        nsamplenotreal=len(self.myglobalparameters.t)

        if apply_window == 1:
            windows = signal.tukey(nsamplenotreal, alpha = 0.05)


        out_dir="temp"
        [algo,swarmsize,maxiter, maxiter_ps]=self.vars_temp_file_5



        # Load fields data
        out_opt_filename = "optim_result"
        out_opt_full_info_filename = f"{out_dir}/{out_opt_filename.split('.')[0]}_full_info.out"

        datacorrection = datalist()


            # =============================================================================
        myvariables = []
        nb_param = len(myvariables)
            
        myVariablesDictionary = {}
        minDict = {}
        maxDict = {}
        totVariablesName = myvariables
            
            
        if fit_delay == 1:
            myVariablesDictionary['delay']=delaymax_guess
            minDict['delay'] = -delay_limit
            maxDict['delay'] =  delay_limit
            totVariablesName = np.append(totVariablesName,'delay')
            
            
        if fit_dilatation == 1:
            tab=[]
            for i in range (0,len(dilatationmax_guess)):                    
                myVariablesDictionary['dilatation '+str(i)]=dilatationmax_guess[i]#leftcoef[count-1]
                minDict['dilatation '+str(i)] = -dilatation_limit[i]
                maxDict['dilatation '+str(i)] = dilatation_limit[i]
                tab = np.append(tab,'dilatation '+str(i))
            totVariablesName = np.append(totVariablesName,tab)
            
            
        if (fit_leftover_noise == 1):
            tab=[]
            for i in range (0,len(leftcoef_guess)):                    
                myVariablesDictionary['leftover '+str(i)]=leftcoef_guess[i]#leftcoef[count-1]
                minDict['leftover '+str(i)] = -leftcoef_limit[i]
                maxDict['leftover '+str(i)] = leftcoef_limit[i]
                tab = np.append(tab,'leftover '+str(i))
            totVariablesName = np.append(totVariablesName,tab)
            ## We take into account the thicknesses and delay as optimization parameters
            # so we put the values and their uncertainty in the corresponding lists
            
            
            #=============================================================================#
            # Instantiate Optimization Problem
            #=============================================================================#


        #*************************************************************************************************************
            # Normalisation
        self.minval = np.array(list(minDict.values()))
        self.maxval = np.array(list(maxDict.values()))
        guess = np.array(list(myVariablesDictionary.values()))

        x0=np.array((guess-self.minval)/(self.maxval-self.minval))
        lb=np.zeros(len(guess))
        up=np.ones(len(guess))

        self.dt=self.myglobalparameters.t.item(2)-self.myglobalparameters.t.item(1)   ## Sample rate

        numberOfTrace = len(data.pulse)
        fopt_init = []
        exposant_ref = 4  # au lieu d'avoir x0 = 0.5 pour la ref, qui est dejà optimal et donc qui fait deconné l'ago d'optim, on aura x0 = 0.5-1e^exposant_ref



        if fit_periodic_sampling: #need an init point for optimization after correction
            # print("Periodic sampling optimization")

            self.mymean = np.mean(data.pulse, axis = 0)

            nu = periodic_sampling_freq_limit*1e12   # 1/s   Hz
            delta_nu = self.myglobalparameters.freq[-1]/(len(self.myglobalparameters.freq)-1) # Hz
            index_nu=int(nu/delta_nu)
            
            maxval_ps = np.array([self.dt/10, 12*2*np.pi*1e12, np.pi])
            minval_ps = np.array([0, 6*2*np.pi*1e12, -np.pi])
            guess_ps = np.array([0,0,0])
            
            x0_ps = (guess_ps-minval_ps)/(maxval_ps-minval_ps)
            lb_ps=np.zeros(len(guess_ps))
            ub_ps=np.ones(len(guess_ps))
            
            def error_periodic(x):
                # x = A, v, phi
                # global mymean
                x = x*(maxval_ps-minval_ps)+minval_ps
                
                ct = x[0]*np.cos(x[1]*self.myglobalparameters.t + x[2])    # s 
                corrected = self.mymean - np.gradient(self.mymean, self.dt)*ct
                
                error = sum(abs((torch_rfft(corrected)[index_nu:])))
                
                #error = 0 # Doesn't work , why?
                #for i in range(index_nu,len(myglobalparameters.freq)):
                  #   error += abs(np.real(np.exp(-j*np.angle(np.fft.rfft(corrected)[i-index_nu])) * np.fft.rfft(corrected)[i]))

                return error
            
            res_ps = optimize.dual_annealing(error_periodic, x0 = x0_ps, maxiter = maxiter_ps, bounds=list(zip(lb_ps, ub_ps)))
            
            xopt_ps = res_ps.x*(maxval_ps-minval_ps)+minval_ps




        if fit_delay or fit_leftover_noise or fit_dilatation:
            # print("Delay and amplitude and dilatation error optimization")
            for trace in range(numberOfTrace) :
                
                # Checking if break button is pressed
                if self.interrupt and self.optimization_process.is_alive():
                    self.optimization_process.terminate()
                    # print("Optimization process terminated")
            
                # print("Time trace "+str(trace))
                self.myinputdata=mydata(data.pulse[trace])    ## We create a variable containing the data related to the measured pulse
                # data.pulse[trace] = [] # why there is this line of code ?
                
                monerreur = self.errorchoice()
                
                # optim = Optimization()
                # objfunc = optim.errorchoice_pyOpt
                
                if fit_leftover_noise:
                    if fit_dilatation:
                        if fit_delay:
                            x0[1] = 0.505 #coef a on evite de commencer l'init à 0 car parfois probleme de convergence
                        else:
                            x0[0] = 0.505  # coef a 
                    elif not fit_dilatation and trace ==0: # si on fit pas la dilatation, on peut utiliser les anciens result d'optim, a part pour la trace 0
                        if fit_delay:
                            x0[1] = 0.505 #coef a on evite de commencer l'init à 0 car parfois probleme de convergence
                        else:
                            x0[0] = 0.505  # coef a 
                
                if trace == reference_number: # on part pas de 0.5 car il diverge vu que c'est la ref elle meme
                    ref_x0= [0.5 - 0.1**exposant_ref]*len(totVariablesName)
                    # on print pas c
                    # if fit_leftover_noise:  
                    #     print('guess')
                    #     print((np.array(ref_x0)*(self.maxval-self.minval)+self.minval)[:-1])
                    #     print('x0')
                    #     print(ref_x0[:-1])
                    # else:
                    #     print('guess')
                    #     print(np.array(ref_x0)*(self.maxval-self.minval)+self.minval)
                    #     print('x0')
                    #     print(ref_x0)
                    # print('errorguess')
                    fopt_init.append(monerreur(ref_x0))
                    # print(fopt_init[-1])
                else:
                    guess= x0*(self.maxval-self.minval)+self.minval
                    # on print seulemnt delay et a , pas c
                    # if fit_leftover_noise:
                    #     print('guess')
                    #     print(guess[:-1])
                    #     print('x0')
                    #     print(x0[:-1])
                    # else:
                    #     print('guess')
                    #     print(guess)
                    #     print('x0')
                    #     print(x0)
                    # print('errorguess')
                    fopt_init.append(monerreur(x0))
                    # print(fopt_init[-1])

                
                
                # ## Optimization dans le cas PyOpt
                # if algo in [1,2,3,4]:
                #     opt_prob = Optimization('Dielectric modeling based on TDS pulse fitting',objfunc)
                #     icount = 0
                #     for nom,varvalue in myVariablesDictionary.items():
                #         #if varvalue>=0:
                #         if trace == reference_number:
                #             opt_prob.addVar(nom,'c',lower = 0,upper = 1,
                #                     value = ref_x0[icount] #normalisation
                #                     )
                #         else:
                #             opt_prob.addVar(nom,'c',lower = 0,upper = 1,
                #                     value = (varvalue-minDict.get(nom))/(maxDict.get(nom)-minDict.get(nom)) #normalisation
                #                     )
                #         icount+=1
                #         #else:
                #         #    opt_prob.addVar(nom,'c',lower = 0,upper = 1,
                #         #                value = -(varvalue-minDict.get(nom))/(maxDict.get(nom)-minDict.get(nom)) #normalisation
                #           #               )    
                #     opt_prob.addObj('f')
                #     #opt_prob.addCon('g1','i') #possibility to add constraints
                #     #opt_prob.addCon('g2','i')
                
                
                # =============================================================================
                # solving the problem with the function in scipy.optimize
                # =============================================================================
                
                
                if  algo==0: 
                    start = time.process_time()
                    xopt,fopt=pso(monerreur,lb,up,swarmsize=swarmsize,minfunc=1e-18,minstep=1e-8,debug=1,phip=0.5,phig=0.5,maxiter=maxiter) ## 'monerreur' function that we want to minimize, 'lb' and 'up' bounds of the problem
                    elapsed_time = time.process_time()-start
                    # print("Time taken by the optimization:",elapsed_time)
                    
                if algo == 5:
                    start = time.process_time()
                    cback=Callback_bfgs(monerreur)
                    if trace == reference_number:
                        res = optimize.minimize(monerreur,ref_x0,method='L-BFGS-B',bounds=list(zip(lb, up)),callback=cback,options={'maxiter':maxiter})
                    else:
                        res = optimize.minimize(monerreur,x0,method='L-BFGS-B',bounds=list(zip(lb, up)),callback=cback,options={'maxiter':maxiter})
                    elapsed_time = time.process_time()-start
                    xopt = res.x
                    fopt = res.fun
                    # print(res.message,"\nTime taken by the optimization:",elapsed_time)
                    
                if algo == 6:
                    start = time.process_time()
                    cback=Callback_slsqp(monerreur)
                    if trace == reference_number:
                        res = optimize.minimize(monerreur,ref_x0,method='SLSQP',bounds=list(zip(lb, up)),callback=cback,options={'maxiter':maxiter, 'ftol': 1e-20})
                    else:
                        res = optimize.minimize(monerreur,x0,method='SLSQP',bounds=list(zip(lb, up)),callback=cback,options={'maxiter':maxiter})
                    elapsed_time = time.process_time()-start
                    xopt = res.x
                    fopt = res.fun
                    # print(res.message,"\nTime taken by the optimization:",elapsed_time)
                    
                if algo==7:
                    start = time.process_time()
                    cback=Callback_annealing(monerreur)
                    res = optimize.dual_annealing(monerreur, bounds=list(zip(lb, up)),callback=cback,maxiter=maxiter)
                    elapsed_time = time.process_time()-start
                    xopt = res.x
                    fopt = res.fun
                    # print(res.message,"\nTime taken by the optimization:",elapsed_time)
                
                
                
                # =============================================================================
                # solving the problem with pyOpt
                # =============================================================================
                
                
                # if  (algo==1)|(algo == 2):
                #     start = time.process_time()
                #     [fopt, xopt, inform] = optimALPSO(opt_prob, swarmsize, maxiter,algo,out_opt_full_info_filename)
                #     elapsed_time = time.process_time()-start
                #     print(inform,"\nTime taken by the optimization:",elapsed_time)
                    
                # if algo ==3:
                #         try:
                #             start = time.process_time()
                #             [fopt, xopt, inform] = optimSLSQP(opt_prob,maxiter)
                #             elapsed_time = time.process_time()-start
                #             print(inform,"\nTime taken by the optimization:",elapsed_time)
                #         except Exception as e:
                #             print(e)
                
                # if algo ==4:
                #         try:
                #             start = time.process_time()
                #             [fopt, xopt, inform] = optimSLSQPpar(opt_prob,maxiter)
                #             elapsed_time = time.process_time()-start
                #             print(inform,"\nTime taken by the optimization:",elapsed_time)
                #         except Exception as e:
                #             print(e)
                  
                # if fit_leftover_noise and not fit_dilatation: 
                #     # si on corrige la dilatation, vaut mieux repartir de 0 sinon divergence
                #     if fit_delay:
                #         x0[1] = xopt[1] #coef a on evite de commencer l'init à 0 car parfois probleme de convergence
                #     else:
                #         x0[0] = xopt[0]  # coef a  
                # =============================================================================
                
                if myrank == 0:
                    xopt = xopt*(self.maxval-self.minval)+self.minval  #denormalize
                    # print(f'The best error was: \t{fopt}')
                    # if(fit_leftover_noise):
                    #     print(f'the best parameters were: \t{xopt[:-1]}\n')
                    # else:
                    #     print(f'the best parameters were: \t{xopt}\n')
                    # =========================================================================
                            
                    myfitteddata=myfitdata(self.myinputdata, xopt)
                    
                    datacorrection.add_trace(myfitteddata.pulse)
                    
                    # =========================================================================
                    # saving the results
                    # ========================================================================
            
                    # result_optimization=[]
                    # result_optimization.append([xopt,fopt]) # instead of result_optimization=[xopt,fopt]
                    result_optimization=[xopt,fopt]
                    if(trace == 0):   # write first time
                        # with open(os.path.join("temp",'temp_file_3.bin'),'wb') as f:
                        #     pickle.dump(result_optimization,f,pickle.HIGHEST_PROTOCOL)
                        result_optimization = [item.tolist() if isinstance(item, np.ndarray) else item for item in result_optimization]
                        result_optimization = [result_optimization]
                        self.vars_temp_file_3 = result_optimization

                    else:  #append after first time
                        # with open(os.path.join("temp",'temp_file_3.bin'),'ab') as f:
                        #     pickle.dump(result_optimization,f,pickle.HIGHEST_PROTOCOL)
                        result_optimization = [item.tolist() if isinstance(item, np.ndarray) else item for item in result_optimization]
                        self.vars_temp_file_3.append(result_optimization) # use .append instead of =



            #TOADD: progressbar   
            ################################### After the optimization loop #############
            
            #Loop for is done : interrupt variable back to False
            self.interrupt=False
            
            if myrank == 0 and not fit_periodic_sampling:  
                datacorrection.moyenne = np.mean(datacorrection.pulse, axis = 0)
                datacorrection.time_std = np.std(datacorrection.pulse, axis = 0)
                datacorrection.freq_std = np.std(torch_rfft(datacorrection.pulse, axis = 1),axis = 0)
                if apply_window == 1:
                    datacorrection.freq_std_with_window = np.std(torch_rfft(datacorrection.pulse*windows, axis = 1),axis = 0)
                #SAVE the result in binary for other modules
                # with open(os.path.join("temp",'temp_file_2.bin'),'wb') as f:
                #     pickle.dump(datacorrection,f,pickle.HIGHEST_PROTOCOL)
                #     pickle.dump(fopt_init,f,pickle.HIGHEST_PROTOCOL)
                self.vars_temp_file_2_datacorrection = datacorrection
                self.vars_temp_file_2_fopt=fopt_init
                
            
                    


        ###################################################
                 #  ****************************************** PERIODIC SAMPLING *******************************************   


        if fit_periodic_sampling:
            # print("Periodic sampling optimization")            

            if fit_delay or fit_leftover_noise or fit_dilatation:
                self.mymean = np.mean(datacorrection.pulse, axis = 0)   
            else:
                self.mymean = np.mean(data.pulse, axis = 0)

            nu = periodic_sampling_freq_limit*1e12   # 1/s   Hz
            delta_nu = self.myglobalparameters.freq[-1]/(len(self.myglobalparameters.freq)-1) # Hz
            index_nu=int(nu/delta_nu)
            
            maxval_ps = np.array([self.dt/10, 12*2*np.pi*1e12, np.pi])
            minval_ps = np.array([0, 6*2*np.pi*1e12, -np.pi])
            guess_ps = xopt_ps
            
            x0_ps = (guess_ps-minval_ps)/(maxval_ps-minval_ps)
            lb_ps=np.zeros(len(guess_ps))
            ub_ps=np.ones(len(guess_ps))
            
            def error_periodic(x):
                # x = A, v, phi
                # global mymean
                x = x*(maxval_ps-minval_ps)+minval_ps
                
                ct = x[0]*np.cos(x[1]*self.myglobalparameters.t + x[2])    # s 
                corrected = self.mymean - np.gradient(self.mymean, self.dt)*ct
                
                error = sum(abs((torch_rfft(corrected)[index_nu:])))
                
                #error = 0
                # for i in range(index_nu,len(myglobalparameters.freq)):
                #     error += abs(np.real(np.exp(-j*np.angle(np.fft.rfft(corrected)[i-index_nu])) * np.fft.rfft(corrected)[i]))

                return error
            
            # print('guess')
            # print(guess_ps)
            # print('x0')
            # print(x0_ps)
            res_ps = optimize.dual_annealing(error_periodic, x0 = x0_ps, maxiter = maxiter_ps, bounds=list(zip(lb_ps, ub_ps)))
            #res_ps = optimize.minimize(error_periodic,x0_ps, method='SLSQP',bounds=list(zip(lb_ps, ub_ps)), options={'maxiter':maxiter_ps})
            #res_ps = optimize.minimize(error_periodic,x0_ps,method='L-BFGS-B',bounds=list(zip(lb_ps, ub_ps)), options={'maxiter':1000})
            #res_ps = pso(error_periodic,lb_ps,ub_ps,swarmsize=100,minfunc=1e-18,minstep=1e-8,debug=1,phip=0.5,phig=0.5,maxiter=100)
            
            xopt_ps = res_ps.x*(maxval_ps-minval_ps)+minval_ps
            fopt_ps = res_ps.fun
            
            result_optimization = [xopt_ps, fopt_ps]
            
            if fit_delay or fit_leftover_noise or fit_dilatation:
                # with open(os.path.join("temp",'temp_file_3.bin'),'ab') as f:
                #     pickle.dump(result_optimization,f,pickle.HIGHEST_PROTOCOL)
                # self.vars_temp_file_3 = result_optimization
                result_optimization = [item.tolist() if isinstance(item, np.ndarray) else item for item in result_optimization]
                self.vars_temp_file_3.append(result_optimization) # use .append instead of =

            else:
                # with open(os.path.join("temp",'temp_file_3.bin'),'wb') as f:
                #     pickle.dump(result_optimization,f,pickle.HIGHEST_PROTOCOL)
                result_optimization = [item.tolist() if isinstance(item, np.ndarray) else item for item in result_optimization]
                self.vars_temp_file_3 = result_optimization
            
            ct = xopt_ps[0]*np.cos(xopt_ps[1]*self.myglobalparameters.t + xopt_ps[2])
            
            if fit_delay or fit_leftover_noise or fit_dilatation:
                for i in range(numberOfTrace):
                    # print("correction of trace {}".format(i))
                    # print(f"correction of trace {i}")
                    datacorrection.pulse[i]= datacorrection.pulse[i] - np.gradient(datacorrection.pulse[i], self.dt)*ct
            else:
                for i in range(numberOfTrace):
                    # print(f"correction of trace {i}")
                    temp = data.pulse[i] - np.gradient(data.pulse[i], self.dt)*ct
                    datacorrection.add_trace(temp)

            
            # print(f'The best error was: \t{fopt_ps}')
            # print(f'the best parameters were: \t{xopt_ps}\n')

            datacorrection.moyenne = np.mean(datacorrection.pulse, axis = 0)
            datacorrection.time_std = np.std(datacorrection.pulse, axis = 0)
            datacorrection.freq_std = np.std(torch_rfft(datacorrection.pulse, axis = 1),axis = 0) 
            if apply_window == 1:
                datacorrection.freq_std_with_window = np.std(torch_rfft(datacorrection.pulse*windows, axis = 1),axis = 0)
                        
            # with open(os.path.join("temp",'temp_file_2.bin'),'wb') as f:
            #     pickle.dump(datacorrection,f,pickle.HIGHEST_PROTOCOL)
            #     pickle.dump(fopt_init,f,pickle.HIGHEST_PROTOCOL)
            self.vars_temp_file_2_datacorrection = datacorrection
            self.vars_temp_file_2_fopt=fopt_init

        ###################################################