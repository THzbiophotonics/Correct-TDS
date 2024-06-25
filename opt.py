# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 12:34:04 2019

@author: nayab, juliettevl
"""
# =============================================================================
# Standard Python modules
# =============================================================================
import os, time
import pickle
from pyswarm import pso   ## Library for optimization
import numpy as np   ## Library to simplify the linear algebra calculations
import scipy.optimize as optimize  ## Library for optimization
import fitf as TDS
import param
from scipy import signal



# =============================================================================
j = 1j


# =============================================================================
# External Python modules (serves for optimization algo #3)
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


#end
# =============================================================================
# Extension modules
# =============================================================================

try:
    from pyOpt import Optimization   ## Library for optimization
    from pyOpt import ALPSO  ## Library for optimization
except:
    if myrank==0:
        print("Error importing pyopt")
    
try:
    from pyOpt import SLSQP  ## Library for optimization
except:
    if myrank==0:
        print("Error importing pyopt SLSQP")


# =============================================================================
# classes we will use
# =============================================================================
        
class myfitdata: 
    def __init__(self, myinput, x):
        self.pulse = fit_input(myinputdata, x)
        self.Spulse = (TDS.torch_rfft(self.pulse))

    

# =============================================================================
class Callback_bfgs(object):
    def __init__(self):
        self.nit = 0
        
    def __call__(self, par, convergence=0):
        self.nit += 1
        with open('algo_bfgs_out.txt', 'a+') as filehandle:
            filehandle.write('\n iteration number %d ; error %s ; parameters %s \r\n' % (self.nit, monerreur(par), par))
            
class Callback_slsqp(object):
    def __init__(self):
        self.nit = 0
        
    def __call__(self, par, convergence=0):
        self.nit += 1
        with open('algo_slsqp_out.txt', 'a+') as filehandle:
            filehandle.write('\n iteration number %d ; error %s ; parameters %s \r\n' % (self.nit, monerreur(par), par))


class Callback_annealing(object):
    def __init__(self):
        self.nit = 0
        
    def __call__(self, par, f, context):
        self.nit += 1
        with open('algo_dualannealing_out.txt', 'a+') as filehandle:
            filehandle.write('\n iteration number %d ; error %s ; parameters %s \r\n' % (self.nit, monerreur(par), par))





#=============================================================================
def errorchoice():
    global myglobalparameters, myinputdata, mode, fit_delay, fit_leftover_noise, fit_dilatation, dt, nsample, maxval, minval
    
    def monerreur(x):
        
        leftover_guess = np.zeros(2)
        delay_guess = 0
        dilatation_coefguess = np.zeros(2)
        
        coef = np.zeros(2) #[a,c]

        x = x*(maxval-minval)+minval
        
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

        Z = np.exp(j*myglobalparameters.w*delay_guess)
        myinputdatacorrected_withdelay = TDS.torch_irfft(Z*myinputdata.Spulse, n = len(myglobalparameters.t))

        leftnoise = np.ones(len(myglobalparameters.t)) - coef[0]*np.ones(len(myglobalparameters.t))   #(1-a)    
        myinputdatacorrected = leftnoise*(myinputdatacorrected_withdelay 
                                          - (dilatation_coefguess[0]*myglobalparameters.t)*np.gradient(myinputdatacorrected_withdelay, dt))
        erreur = np.linalg.norm(myreferencedata.Pulseinit[:nsample] - myinputdatacorrected[:nsample] )/np.linalg.norm(myreferencedata.Pulseinit[:nsample])
        return erreur
            
    return monerreur


def fit_input(myinputdata, x):
    global myglobalparameters, mode, fit_delay, fit_leftover_noise, fit_dilatation, dt
    
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

    Z = np.exp(j*myglobalparameters.w*delay_guess)
    myinputdatacorrected_withdelay = TDS.torch_irfft(Z*myinputdata.Spulse, n = len(myglobalparameters.t))

    leftnoise = np.ones(len(myglobalparameters.t)) - coef[0]*np.ones(len(myglobalparameters.t))   #(1-a)    
    myinputdatacorrected = leftnoise*(myinputdatacorrected_withdelay  
                                      - (dilatation_coefguess[0]*myglobalparameters.t)*np.gradient(myinputdatacorrected_withdelay, dt))

    return myinputdatacorrected


# =============================================================================
def errorchoice_pyOpt(): 
    def objfunc(x):  ## Function used in the Optimization function from pyOpt. For more details see http://www.pyopt.org/quickguide/quickguide.html
        monerreur = errorchoice()
        f = monerreur(x)
        fail = 0
        return f, 1, fail
    return objfunc


# =============================================================================

def optimALPSO(opt_prob, swarmsize, maxiter,algo,out_opt_full_info_filename): #ok
    if algo == 2:
        alpso_none = ALPSO(pll_type='SPM')
    else:
        alpso_none = ALPSO()
    #alpso_none.setOption('fileout',1)
    #alpso_none.setOption('filename',out_opt_full_info_filename)
    alpso_none.setOption('SwarmSize',swarmsize)
    alpso_none.setOption('maxInnerIter',6)
    alpso_none.setOption('etol',1e-5)
    alpso_none.setOption('rtol',1e-10)
    alpso_none.setOption('atol',1e-10)
    alpso_none.setOption('vcrazy',1e-4)
    alpso_none.setOption('dt',1e0)
    alpso_none.setOption('maxOuterIter',maxiter)
    alpso_none.setOption('stopCriteria',0)#Stopping Criteria Flag (0 - maxIters, 1 - convergence)
    alpso_none.setOption('printInnerIters',1)               
    alpso_none.setOption('printOuterIters',1)
    alpso_none.setOption('HoodSize',int(swarmsize/100))
    return(alpso_none(opt_prob))
    
def optimSLSQP(opt_prob,maxiter):#ok
    slsqp_none = SLSQP()
    #slsqp_none.setOption('IPRINT',1)
    #slsqp_none.setOption('IFILE',out_opt_full_info_filename)
    slsqp_none.setOption('MAXIT',maxiter)
    slsqp_none.setOption('IOUT',15) 
    slsqp_none.setOption('ACC',1e-20)
    return(slsqp_none(opt_prob))

def optimSLSQPpar(opt_prob,maxiter): # arecopierdansdoublet #ok
          
    slsqp_none = SLSQP() # arecopierdansdoublet

    #slsqp_none.setOption('IPRINT',1)
    #slsqp_none.setOption('IFILE',out_opt_full_info_filename)
    slsqp_none.setOption('MAXIT',maxiter)
    slsqp_none.setOption('IOUT',12) 
    slsqp_none.setOption('ACC',1e-24)
    return(slsqp_none(opt_prob,sens_mode='pgc')) 

            

# =============================================================================
# We load the model choices
# =============================================================================
f=open(os.path.join("temp",'temp_file_1_ini.bin'),'rb') #ok
[path_data, path_data_ref, reference_number, fit_dilatation, dilatation_limit, dilatationmax_guess, 
 freqWindow, timeWindow, fit_delay, delaymax_guess, delay_limit, mode, nsample,
 fit_periodic_sampling, periodic_sampling_freq_limit, fit_leftover_noise, leftcoef_guess, leftcoef_limit]=pickle.load(f)
f.close()

data = TDS.datalist()
f=open(os.path.join("temp",'temp_file_6.bin'),'rb')
data = pickle.load(f)
myreferencedata=pickle.load(f) # champs
f.close()

myglobalparameters = TDS.globalparameters()
f=open(os.path.join("temp",'temp_file_7.bin'),'rb')
myglobalparameters = pickle.load(f)
apply_window = pickle.load(f)
nsamplenotreal=len(myglobalparameters.t)
f.close()

if apply_window == 1:
    windows = param.window(nsamplenotreal)


out_dir="temp"

f=open(os.path.join("temp",'temp_file_5.bin'),'rb')
[algo,swarmsize,maxiter, maxiter_ps]=pickle.load(f)
f.close()



# Load fields data
out_opt_filename = "optim_result"
out_opt_full_info_filename=os.path.join(out_dir,'{0}_full_info.out'.format(out_opt_filename.split('.')[0]))

datacorrection = TDS.datalist()


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
minval = np.array(list(minDict.values()))
maxval = np.array(list(maxDict.values()))
guess = np.array(list(myVariablesDictionary.values()))

x0=np.array((guess-minval)/(maxval-minval))
lb=np.zeros(len(guess))
up=np.ones(len(guess))

dt=myglobalparameters.t.item(2)-myglobalparameters.t.item(1)   ## Sample rate

numberOfTrace = len(data.pulse)
fopt_init = []
exposant_ref = 4  # au lieu d'avoir x0 = 0.5 pour la ref, qui est dejà optimal et donc qui fait deconné l'ago d'optim, on aura x0 = 0.5-1e^exposant_ref



if fit_periodic_sampling: #need an init point for optimization after correction
    print("Periodic sampling optimization")

    mymean = np.mean(data.pulse, axis = 0)

    nu = periodic_sampling_freq_limit*1e12   # 1/s   Hz
    delta_nu = myglobalparameters.freq[-1]/(len(myglobalparameters.freq)-1) # Hz
    index_nu=int(nu/delta_nu)
    
    maxval_ps = np.array([dt/10, 12*2*np.pi*1e12, np.pi])
    minval_ps = np.array([0, 6*2*np.pi*1e12, -np.pi])
    guess_ps = np.array([0,0,0])
    
    x0_ps = (guess_ps-minval_ps)/(maxval_ps-minval_ps)
    lb_ps=np.zeros(len(guess_ps))
    ub_ps=np.ones(len(guess_ps))
    
    def error_periodic(x):
        # x = A, v, phi
        global mymean
        x = x*(maxval_ps-minval_ps)+minval_ps
        
        ct = x[0]*np.cos(x[1]*myglobalparameters.t + x[2])    # s 
        corrected = mymean - np.gradient(mymean, dt)*ct
        
        error = sum(abs((TDS.torch_rfft(corrected)[index_nu:])))
        
        #error = 0 # Doesn't work , why?
        #for i in range(index_nu,len(myglobalparameters.freq)):
         #   error += abs(np.real(np.exp(-j*np.angle(np.fft.rfft(corrected)[i-index_nu])) * np.fft.rfft(corrected)[i]))

        return error
    
    res_ps = optimize.dual_annealing(error_periodic, x0 = x0_ps, maxiter = maxiter_ps, bounds=list(zip(lb_ps, ub_ps)))
    
    xopt_ps = res_ps.x*(maxval_ps-minval_ps)+minval_ps




if fit_delay or fit_leftover_noise or fit_dilatation:
    print("Delay and amplitude and dilatation error optimization")
    for trace in  range(numberOfTrace) :
        print("Time trace "+str(trace))
        myinputdata=TDS.mydata(data.pulse[trace])    ## We create a variable containing the data related to the measured pulse
        data.pulse[trace] = []
        
        monerreur = errorchoice()
        objfunc = errorchoice_pyOpt()
        
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
            if fit_leftover_noise:  
                print('guess')
                print((np.array(ref_x0)*(maxval-minval)+minval)[:-1])
                print('x0')
                print(ref_x0[:-1])
            else:
                print('guess')
                print(np.array(ref_x0)*(maxval-minval)+minval)
                print('x0')
                print(ref_x0)
            print('errorguess')
            fopt_init.append(monerreur(ref_x0))
            print(fopt_init[-1])
        else:
            guess= x0*(maxval-minval)+minval
            # on print seulemnt delay et a , pas c
            if fit_leftover_noise:
                print('guess')
                print(guess[:-1])
                print('x0')
                print(x0[:-1])
            else:
                print('guess')
                print(guess)
                print('x0')
                print(x0)
            print('errorguess')
            fopt_init.append(monerreur(x0))
            print(fopt_init[-1])

        
        
        ## Optimization dans le cas PyOpt
        if algo in [1,2,3,4]:
            opt_prob = Optimization('Dielectric modeling based on TDS pulse fitting',objfunc)
            icount = 0
            for nom,varvalue in myVariablesDictionary.items():
                #if varvalue>=0:
                if trace == reference_number:
                    opt_prob.addVar(nom,'c',lower = 0,upper = 1,
                            value = ref_x0[icount] #normalisation
                            )
                else:
                    opt_prob.addVar(nom,'c',lower = 0,upper = 1,
                            value = (varvalue-minDict.get(nom))/(maxDict.get(nom)-minDict.get(nom)) #normalisation
                            )
                icount+=1
                #else:
                #    opt_prob.addVar(nom,'c',lower = 0,upper = 1,
                #                value = -(varvalue-minDict.get(nom))/(maxDict.get(nom)-minDict.get(nom)) #normalisation
                 #               )    
            opt_prob.addObj('f')
            #opt_prob.addCon('g1','i') #possibility to add constraints
            #opt_prob.addCon('g2','i')
        
        
        # =============================================================================
        # solving the problem with the function in scipy.optimize
        # =============================================================================
        
        
        if  algo==0: 
            start = time.process_time()
            xopt,fopt=pso(monerreur,lb,up,swarmsize=swarmsize,minfunc=1e-18,minstep=1e-8,debug=1,phip=0.5,phig=0.5,maxiter=maxiter) ## 'monerreur' function that we want to minimize, 'lb' and 'up' bounds of the problem
            elapsed_time = time.process_time()-start
            print("Time taken by the optimization:",elapsed_time)
            
        if algo == 5:
            start = time.process_time()
            cback=Callback_bfgs()
            if trace == reference_number:
                res = optimize.minimize(monerreur,ref_x0,method='L-BFGS-B',bounds=list(zip(lb, up)),callback=cback,options={'maxiter':maxiter})
            else:
                res = optimize.minimize(monerreur,x0,method='L-BFGS-B',bounds=list(zip(lb, up)),callback=cback,options={'maxiter':maxiter})
            elapsed_time = time.process_time()-start
            xopt = res.x
            fopt = res.fun
            print(res.message,"\nTime taken by the optimization:",elapsed_time)
            
        if algo == 6:
            start = time.process_time()
            cback=Callback_slsqp()
            if trace == reference_number:
                res = optimize.minimize(monerreur,ref_x0,method='SLSQP',bounds=list(zip(lb, up)),callback=cback,options={'maxiter':maxiter, 'ftol': 1e-20})
            else:
                res = optimize.minimize(monerreur,x0,method='SLSQP',bounds=list(zip(lb, up)),callback=cback,options={'maxiter':maxiter})
            elapsed_time = time.process_time()-start
            xopt = res.x
            fopt = res.fun
            print(res.message,"\nTime taken by the optimization:",elapsed_time)
            
        if algo==7:
            start = time.process_time()
            cback=Callback_annealing()
            res = optimize.dual_annealing(monerreur, bounds=list(zip(lb, up)),callback=cback,maxiter=maxiter)
            elapsed_time = time.process_time()-start
            xopt = res.x
            fopt = res.fun
            print(res.message,"\nTime taken by the optimization:",elapsed_time)
        
        
        
        # =============================================================================
        # solving the problem with pyOpt
        # =============================================================================
        
        
        if  (algo==1)|(algo == 2):
            start = time.process_time()
            [fopt, xopt, inform] = optimALPSO(opt_prob, swarmsize, maxiter,algo,out_opt_full_info_filename)
            elapsed_time = time.process_time()-start
            print(inform,"\nTime taken by the optimization:",elapsed_time)
            
        if algo ==3:
                try:
                    start = time.process_time()
                    [fopt, xopt, inform] = optimSLSQP(opt_prob,maxiter)
                    elapsed_time = time.process_time()-start
                    print(inform,"\nTime taken by the optimization:",elapsed_time)
                except Exception as e:
                    print(e)
        
        if algo ==4:
                try:
                    start = time.process_time()
                    [fopt, xopt, inform] = optimSLSQPpar(opt_prob,maxiter)
                    elapsed_time = time.process_time()-start
                    print(inform,"\nTime taken by the optimization:",elapsed_time)
                except Exception as e:
                    print(e)
          
        if fit_leftover_noise and not fit_dilatation: 
            # si on corrige la dilatation, vaut mieux repartir de 0 sinon divergence
            if fit_delay:
                x0[1] = xopt[1] #coef a on evite de commencer l'init à 0 car parfois probleme de convergence
            else:
                x0[0] = xopt[0]  # coef a  
        # =============================================================================
        
        if myrank == 0:
            xopt = xopt*(maxval-minval)+minval  #denormalize
            print('The best error was: \t{}'.format(fopt))
            if(fit_leftover_noise):
                print('the best parameters were: \t{}\n'.format(xopt[:-1]))
            else:
                print('the best parameters were: \t{}\n'.format(xopt))
            # =========================================================================
                    
            myfitteddata=myfitdata(myinputdata, xopt)
            
            datacorrection.add_trace(myfitteddata.pulse)
            # =========================================================================
            # saving the results
            # ========================================================================
    
            result_optimization=[xopt,fopt]
            if(trace == 0):   # write first time
                f=open(os.path.join("temp",'temp_file_3.bin'),'wb')
                pickle.dump(result_optimization,f,pickle.HIGHEST_PROTOCOL)
                f.close()
            else:  #append after first time
                f=open(os.path.join("temp",'temp_file_3.bin'),'ab')
                pickle.dump(result_optimization,f,pickle.HIGHEST_PROTOCOL)
                f.close()
                   
    ################################### After the optimization loop #############
    
    
    if myrank == 0 and not fit_periodic_sampling:  
        datacorrection.moyenne = np.mean(datacorrection.pulse, axis = 0)
        datacorrection.time_std = np.std(datacorrection.pulse, axis = 0)
        datacorrection.freq_std = np.std(TDS.torch_rfft(datacorrection.pulse, axis = 1),axis = 0)
        if apply_window == 1:
            datacorrection.freq_std_with_window = np.std(TDS.torch_rfft(datacorrection.pulse*windows, axis = 1),axis = 0)
        #SAVE the result in binary for other modules
        f=open(os.path.join("temp",'temp_file_2.bin'),'wb')
        pickle.dump(datacorrection,f,pickle.HIGHEST_PROTOCOL)
        pickle.dump(fopt_init,f,pickle.HIGHEST_PROTOCOL)
        f.close()
        
    
            


###################################################
         #  ****************************************** PERIODIC SAMPLING *******************************************   


if fit_periodic_sampling:
    print("Periodic sampling optimization")

    if fit_delay or fit_leftover_noise or fit_dilatation:
        mymean = np.mean(datacorrection.pulse, axis = 0)   
    else:
        mymean = np.mean(data.pulse, axis = 0)
    if param.apply_window_PS:
        mymean *= param.window_PS(len(mymean))
    nu = periodic_sampling_freq_limit*1e12   # 1/s   Hz
    delta_nu = myglobalparameters.freq[-1]/(len(myglobalparameters.freq)-1) # Hz
    index_nu=int(nu/delta_nu)
    
    maxval_ps = np.array([dt/10, 12*2*np.pi*1e12, np.pi])
    minval_ps = np.array([0, 6*2*np.pi*1e12, -np.pi])
    guess_ps = xopt_ps
    
    x0_ps = (guess_ps-minval_ps)/(maxval_ps-minval_ps)
    lb_ps=np.zeros(len(guess_ps))
    ub_ps=np.ones(len(guess_ps))
    
    def error_periodic(x):
        # x = A, v, phi
        global mymean
        x = x*(maxval_ps-minval_ps)+minval_ps
        
        ct = x[0]*np.cos(x[1]*myglobalparameters.t + x[2])    # s 
        corrected = mymean - np.gradient(mymean, dt)*ct
        
        error = sum(abs((TDS.torch_rfft(corrected)[index_nu:])))
        
        #error = 0
        #for i in range(index_nu,len(myglobalparameters.freq)):
         #   error += abs(np.real(np.exp(-j*np.angle(np.fft.rfft(corrected)[i-index_nu])) * np.fft.rfft(corrected)[i]))

        return error
    
    print('guess')
    print(guess_ps)
    print('x0')
    print(x0_ps)
    res_ps = optimize.dual_annealing(error_periodic, x0 = x0_ps, maxiter = maxiter_ps, bounds=list(zip(lb_ps, ub_ps)))
    #res_ps = optimize.minimize(error_periodic,x0_ps, method='SLSQP',bounds=list(zip(lb_ps, ub_ps)), options={'maxiter':maxiter_ps})
    #res_ps = optimize.minimize(error_periodic,x0_ps,method='L-BFGS-B',bounds=list(zip(lb_ps, ub_ps)), options={'maxiter':1000})
    #res_ps = pso(error_periodic,lb_ps,ub_ps,swarmsize=100,minfunc=1e-18,minstep=1e-8,debug=1,phip=0.5,phig=0.5,maxiter=100)
    
    xopt_ps = res_ps.x*(maxval_ps-minval_ps)+minval_ps
    fopt_ps = res_ps.fun
    
    result_optimization = [xopt_ps, fopt_ps]
    
    if fit_delay or fit_leftover_noise or fit_dilatation:
        f=open(os.path.join("temp",'temp_file_3.bin'),'ab')
    else:
        f=open(os.path.join("temp",'temp_file_3.bin'),'wb')
    pickle.dump(result_optimization,f,pickle.HIGHEST_PROTOCOL)
    f.close()
    
    ct = xopt_ps[0]*np.cos(xopt_ps[1]*myglobalparameters.t + xopt_ps[2])
    
    if fit_delay or fit_leftover_noise or fit_dilatation:
        for i in range(numberOfTrace):
            print("correction of trace {}".format(i))
            datacorrection.pulse[i]= datacorrection.pulse[i] - np.gradient(datacorrection.pulse[i], dt)*ct
    else:
        for i in range(numberOfTrace):
            print("correction of trace {}".format(i))
            temp = data.pulse[i] - np.gradient(data.pulse[i], dt)*ct
            datacorrection.add_trace(temp)

    
    print('The best error was: \t{}'.format(fopt_ps))
    print('the best parameters were: \t{}\n'.format(xopt_ps))

    datacorrection.moyenne = np.mean(datacorrection.pulse, axis = 0)
    datacorrection.time_std = np.std(datacorrection.pulse, axis = 0)
    datacorrection.freq_std = np.std(TDS.torch_rfft(datacorrection.pulse, axis = 1),axis = 0) 
    if apply_window == 1:
        datacorrection.freq_std_with_window = np.std(TDS.torch_rfft(datacorrection.pulse*windows, axis = 1),axis = 0)
                
    f=open(os.path.join("temp",'temp_file_2.bin'),'wb')
    pickle.dump(datacorrection,f,pickle.HIGHEST_PROTOCOL)
    pickle.dump(fopt_init,f,pickle.HIGHEST_PROTOCOL)
    f.close()

###################################################


