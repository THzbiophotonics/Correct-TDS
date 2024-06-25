from scipy import signal

apply_window = 1                                                        #For display and saved transfer function
window = lambda nsample: signal.windows.tukey(nsample, alpha = 0.05)    #window is a FUNCTION taking the number of sample point as a parameter

apply_window_PS = 1                                                     #For periodic sampling correction
window_PS = window