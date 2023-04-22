# encoding: utf-8
# main run_script for mode lock program
# sourcery skip: bin-op-identity, hoist-statement-from-loop
# the dircetory of the program is d:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/data/
'''
created by ksr in 2023/3/21
this script is used to simulate the mode lock system
include the following parts:
1. the parameters of the system
2. the optical loop
3. the optoelectronic loop
4. the output image
'''
import csv
import math
import time as time1

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy as sp
from tqdm import tqdm, trange

from components_simplify.Electrical_components import *
from components_simplify.Optical_components import *
from methods_simplify import *
from methods_simplify.Electrical_functions import *
from methods_simplify.IP_CQEM_FD import IP_CQEM_FD
from methods_simplify.Optical_functions import *

' the parameters of the system '
# region the parameters of the system

n_re = 1.45 # the refractive index of the cladding
c = 299792.458 # the speed of light nm/ps
nt = 2**12 # the number of the time samples
f_n = 10e9 # the need frequency 10GHz
wave_length = 1550 # the wave length of the laser 1550nm
time = 120 # ps
twindow = time 
dt = time/nt
t = np.linspace(-twindow/2, twindow/2-dt, nt) # the time sequence
df = 1/(nt*dt) # the frequency step
f = np.linspace(-nt/2, nt/2-1, nt)*df # the frequency sequence
w = 2*np.pi*f # the angular frequency sequence
f_w = 2*np.pi*c/wave_length # the frequency of the laser
lambda_1 = c/(f+c/wave_length) # the wavelength sequence
fo = c/wave_length # the center frequency of the laser
dz = 0.00001 # longitudinal step (km)
tol = 2e-4 # tolerance for the convergence
P0 = 0 # the power of the laser
u0_1 = math.sqrt(P0)*np.ones(nt)
u0 = u0_1+Aoppm(t, 2*np.pi*fo) # the initial signal
# import the data from matlab , using it as comparition
# mat_loc = 'D:/Onedrive/OneDrive - CQU/Code/Python/Activate_Mode_Lock/data/u0_matlab.mat'
# load_mat = sp.io.loadmat(mat_loc)
# u0 = load_mat['u0']
# u0 = np.reshape(u0, (4096,))

# fiber_prameters = [id, length, diameter, aeff, gamma, alpha, beta2, n2]
# SM_1550
SM_1550 = fiber('SMF', 0.010, 1, 60, 1.5539, 0.2/4.343, -22, 2.3)
# SM_1550.gamma = 2*np.pi*SM_1550.n2/wave_length/SM_1550.aeff*1e4 # the gamma of the fiber
SM_1550.gamma = 1.4
SM_1550.alpha = 0/4.343
SM_1550.raman = 0 #exculde the raman effect
SM_1550.ssp = 0 # disable the self phase modulation
SM_1550.betaw = [0, 0, -22, 0.126]

# EDF_1550
EDF_1550 = fiber('EDF', 0.0001, 1, 60, 3.69, 0.2/4.343, -0.130, 2.3)
EDF_1550.length = 0.0001
EDF_1550.gamma = 2*np.pi*SM_1550.n2/wave_length/SM_1550.aeff*1e4 # the gamma of the fiber
EDF_1550.raman = 0 #exculde the raman effect
EDF_1550.ssp = 0 # disable the self phase modulation
EDF_1550.gssdB = 46 # the gain
EDF_1550.EsatpJ = 10 # the saturation PJ
EDF_1550.PsatdBm = 10*math.log10(EDF_1550.EsatpJ/time/0.001) # the saturation power dBm
EDF_1550.lambda_gain = wave_length
EDF_1550.lambda_bw = 40 # the bandwidth of the gain
EDF_1550.fc = c/EDF_1550.lambda_gain # the center frequency of gain
EDF_1550.fbw = c/(EDF_1550.lambda_gain)**2*EDF_1550.lambda_bw # the bandwidth of gain
EDF_1550.betaw = [0, 0, -22.104, 0.171]


#DCF_1550
DCF_1550 = fiber('DCF', 0.0005, 1, 60, 3.2, 0, 123.523, 2.3)
DCF_1550.betaw = [0, 0, 123.523, 0]

# DSF_1550
DSF_1550 = fiber('DSF', 0.090, 1, 60, 3.84, 0, 0, 2.3)
DSF_1550.betaw = [0, 0, 3.523, -0.374]

# SA modulator parameters = [name, q_0, q_1, P_0, Peak_Power]
SA_1550 = Absorber('SA', 0.05, 0.1, 1e2, max(u0**2))
# SA_1550 = Absorber('SA', 0.05, 0.1, 1e2, max(u0**2))

# optical oc parameters = [rohout]
opti_oc = OptiOc(0.01) # the optical circulator

# optical filter parameters = [lamda_c, lamda_bw, n]
opti_filter = OpticalFilter(1550, 10, 1) # the optical filter

# modulator parameters
Modulator_P = Modulator('Phase_Modulator', 3.5, 0, 180/180*math.pi)
Modulator_A = Modulator('Amplitude_Modulator', 6, 3, 180/180*math.pi)

# # input sin siganl
f_Hz = 100e6
f_Thz = f_Hz/1e12
Vp = 10 # Modultor siganl Vp
# M_Signal_in = sin_signal_generator(1, f_Hz, t, (0/180)*math.pi )

# endregion


' the optical loop '
# region the optical loop
# round 500 trips
u = u0
N_trip = 500
spec_z = np.empty((N_trip, nt))
u_z = np.empty((N_trip, nt))
watch_parameters = np.array([])
print('the optical loop is running...')
for _ in trange(N_trip):
    # time1.sleep(0.01)
    u, nf, Plotdata_EDF = IP_CQEM_FD(u,dt,dz,EDF_1550,fo,tol,1,1) # across the EDF
    # u, nf, Plotdata_SM = IP_CQEM_FD(u,dt,dz,SM_1550,fo,tol,1,1)  # across the SMF
    # u, nf, Plotdata_DCF = IP_CQEM_FD(u,dt,dz,DCF_1550,fo,tol,1,1) # across the DCF
    # u, nf, Plotdata_DSF = IP_CQEM_FD(u,dt,dz,DSF_1550,fo,tol,1,1) # across the DSF

    # u = u*Saturation(u, SA_1550.q_0, SA_1550.q_1, SA_1550.P_0) # the SA modulation, waiting change it into PM modulation
    # u = u*Modulation_way_1(f_Thz, Modulator_A.Gamb, Modulator_A.Gamm, t, Modulator_A.phase)
    # u = u*Modulation_way_2(f_Thz, Modulator_P.phase, t, Modulator_P.Vpi, Vp)
    u = Modulation_way_simple(u, 3.5, 3.5, 3.5, nt, 10, 7, t)
    uout, ur = opti_coupler(u, 0, opti_oc.rho_out)
    u = ur
    u_OC = u
    u = filter_gauss(u, opti_filter.f3dB, opti_filter.fc, opti_filter.n, fo, df)
    
    spec = abs(Et2Ef(uout, twindow)**2)
    specnorm = spec/lambda_1**2
    spec_z[_:] = specnorm
    u_z[_:] = np.real(uout)

# u = u*Modulation_way_1(f_Thz, Modulator_A.Gamb, Modulator_A.Gamm, t, Modulator_A.phase)
# u = u*Modulation_way_2(f_Thz, Modulator_P.phase, t, Modulator_P.Vpi, Vp)

# endregion
# watch_parameters = np.append(watch_parameters, Saturation(u, SA_1550.q_0, SA_1550.q_1, SA_1550.P_0), axis=0)
# with open('D:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/z_data/watch_parameters.txt', 'w') as files_1:
#     write_data = csv.writer(files_1)
#     write_data.writerow(watch_parameters)
# region the plot ways
# plt.ion()
plt.figure(1)
plt.plot(t, np.abs(uout)**2, label='uout', color='g')
plt.xlabel('Time/ps')
plt.title('The Output Signal')
plt.ylabel('Amplifier/au')
plt.grid(axis='both', linewidth=0.5)
plt.savefig('D:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/z_data/Time_domain_output.png')

plt.figure(2)
plt.plot(c/(f+fo), 10*np.log10(specnorm), label='uout')
plt.grid(axis='both', linewidth=0.5)
plt.xlabel('Wavelength/nm')
plt.ylabel('Spectrum/dB')
plt.title('The Output Spectrum')
# plt.xlim(1500, 1600)
# plt.ylim(-180, 0)
plt.savefig('D:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/z_data/Spectrum_output.png')

# fig3, ax3 = plt.subplots(1, 2)
# fig3.set_size_inches(10, 5.5)
# c1 = ax3[0].pcolor(t, np.arange(N_trip), np.abs(u_z)**2, cmap='jet')
# ax3[0].set_title('The Output Signal Field Evolution')
# ax3[0].set_xlabel('Time/ps')
# ax3[0].set_ylabel('Trips')
# fig3.colorbar(c1, ax=ax3[0])
# c1 = ax3[1].pcolor(c/(f+fo), np.arange(N_trip), spec_z, cmap='jet')
# ax3[1].set_title('The Output Spectrum Field Evolution')
# ax3[1].set_xlabel('Wavelength/nm')
# ax3[1].set_ylabel('Trips')
# fig3.colorbar(c1, ax=ax3[1])

# fig4, ax4 = plt.subplots(1, 2)
# fig4.set_size_inches(10, 5.5)
# spec = Et2Ef(uout, twindow)**2
# specnorm = spec/lambda_1**2
# c = ax4[0].plot(c/(f+fo), 10*np.log10(specnorm), label='uout')
# ax4[0].set_xlabel('Wavelength/nm')
# ax4[0].set_ylabel('Spectrum/dB')
# ax4[0].set_title('The Output Spectrum')

# c = ax4[1].plot(t, np.abs(uout)**2, color = 'g')
# c = ax4[1].plot(t, np.abs(u0)**2, color = 'b')
# ax4[1].set_xlabel('Time/ps')
# ax4[1].set_ylabel('Peak Power/W')
# ax4[1].set_title('Initial (blue) and Final (green) Pulse Shapes')

# # chirp calculation
# phase_out = np.angle(uout)
# delta_w = np.diff(phase_out)
# Eout = np.abs(uout)**2
# uout_max, uout_max_index = max(Eout), np.argmax(Eout)
# width, I_1, I_2 = fwhm(Eout)
# w_range = np.arange(I_1, I_2, 1)
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(t, Eout, 'g-')
# ax1.set_ylabel('Eout')
# ax1.set_xlabel('Time/ps')
# ax1.legend(loc='best')

# ax2.plot(t[w_range], delta_w[w_range], 'b-')
# ax2.set_ylabel('Eout_half')
# ax2.legend(loc='best')
# ax2.set_title('The Chirp')

# endregion



' the optoelectronic loop '
# region the optoelectronic loop





# endregion



' the output result '
# region the output result

basic_freq = 3e8/n_re/((SM_1550.length+EDF_1550.length+DCF_1550.length)+DSF_1550.length)
all_beta2 = SM_1550.length*SM_1550.beta2 + EDF_1550.length*EDF_1550.beta2 + DCF_1550.length*DCF_1550.beta2 + DSF_1550.length*DSF_1550.beta2
all_l = (SM_1550.length + EDF_1550.length + DCF_1550.length + DSF_1550.length)*1000
print('output peak power is: ', np.max(np.abs(uout)**2), 'W')
print('output pulse energy in pJ is:', dt*(np.sum(np.abs(uout)**2)), 'pJ')
print('the all length is: %.2f m' % all_l)
print('the all beta2 is: ', all_beta2, 'ps^2')
print('the basic frequency is: ', basic_freq/1e9 , 'MHz')
# print(f'the output result is: {basic_freq}')

plt.show()


# endregion