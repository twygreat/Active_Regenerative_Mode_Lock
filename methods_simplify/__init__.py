# basic instruction for the mode lock
# include some basic functions

import numpy as np
import matplotlib.pyplot as plt


# power changes functions

def power_dBm_to_W(power_dBm):
    return 10**(power_dBm/10.0)*1e-3

def power_W_to_dBm(power_W):
    return 10*np.log10(power_W*1e3)

def power_dBm_to_mW(power_dBm):
    return 10**(power_dBm/10.0)

def power_mW_to_dBm(power_mW):
    return 10*np.log10(power_mW)

# some assistant functions
def sin_signal_generator(A, f, t, phi):
    return A*np.sin(2*np.pi*f*t+phi)

