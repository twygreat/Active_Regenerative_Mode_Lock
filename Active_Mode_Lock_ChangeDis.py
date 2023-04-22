import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fft
import math
from tqdm import trange
from methods_simplify.IP_CQEM_FD import IP_CQEM_FD
from methods_simplify.Optical_functions import *
from components_simplify.Optical_components import *
import pandas as pd

""" This code simplfilied descirbe the pulse of active mode lock laser """

# Basic parameters
tl = -300
tr = 300
T = tr - tl
N = 2**14
dt = T/N
t = np.linspace(tl, tr-dt, N)
f = np.linspace(-N/2, N/2-1, N)*1/T
w = 2*np.pi*f
c = 299792.458 # ps/nm
wavelength = 1550
fo = c/wavelength
lambda_1 = c/(f+c/wavelength)
freq_all = f + fo
lambda_all = c/freq_all
omiga_all = 2*np.pi/T*np.linspace(-N/2, N/2-1, N)
df = 1/(N*dt)

# IM parameters
mf = 10
T_AM = 1000/mf
b = 0
m = 1
w_intense = 2*np.pi*mf/1000
nn = min(2*round(500/mf*(N/T)), N)
piecewise_1 = np.zeros(int((N-nn)/2))
piecewise_2 = np.ones(int(nn))
piecewise_3 = np.zeros(int((N-nn)/2))
piecewise = np.concatenate((piecewise_1, piecewise_2, piecewise_3))
T_intense = 0.5*(1+np.sin(np.pi/2*(b+m*np.cos(w_intense*t))))*piecewise
b_phase = 0
m_phase = 10
w_phase = 2*np.pi*mf/1000
T_phase = np.exp(1j*(b_phase+m_phase*np.cos(w_phase*t)))*piecewise

# optical oc parameters = [rohout]
opti_oc = OptiOc(0.01) # the optical circulator

# optical filter parameters = [lamda_c, lamda_bw, n]
opti_filter = OpticalFilter(1550, 10, 1) # the optical filter

# EDF
EDF_1550 = fiber('EDF', 0.0001, 1, 60, 3.69, 0.2/4.343, -0.130, 2.3)
EDF_1550.length = 0.0005
EDF_1550.gamma = 2*np.pi*EDF_1550.n2/wavelength/EDF_1550.aeff*1e4 # the gamma of the fiber
EDF_1550.raman = 0 #exculde the raman effect
EDF_1550.ssp = 0 # disable the self phase modulation
EDF_1550.gssdB = 46 # the gain
EDF_1550.EsatpJ = 10 # the saturation PJ
EDF_1550.PsatdBm = 10*math.log10(EDF_1550.EsatpJ/T/0.001) # the saturation power dBm
EDF_1550.lambda_gain = wavelength
EDF_1550.lambda_bw = 40 # the bandwidth of the gain
EDF_1550.fc = c/EDF_1550.lambda_gain # the center frequency of gain
EDF_1550.fbw = c/(EDF_1550.lambda_gain)**2*EDF_1550.lambda_bw # the bandwidth of gain
EDF_1550.betaw = [0, 0, 28.040448, 0.171]

# SM_1550
SM_1550 = fiber('SMF', 0.010, 1, 60, 1.5539, 0.2/4.343, -22, 2.3)
# SM_1550.gamma = 2*np.pi*SM_1550.n2/wave_length/SM_1550.aeff*1e4 # the gamma of the fiber
SM_1550.gamma = 1.4
SM_1550.alpha = 0/4.343
SM_1550.raman = 0 #exculde the raman effect
SM_1550.ssp = 0 # disable the self phase modulation
SM_1550.betaw = [0, 0, -22, 0.126]

#DCF_1550
DCF_1550 = fiber('DCF', 0.0005, 1, 60, 3.2, 0, 123.523, 2.3)
DCF_1550.betaw = [0, 0, 123.523, 0]

# DSF_1550
DSF_1550 = fiber('DSF', 0.090, 1, 60, 3.84, 0, 0, 2.3)
DSF_1550.betaw = [0, 0, 3.523, -0.374]

g0 = 2
Psat = 50/1000
L1 = 0.015
L2 = 0.0001
L3 = 0.035
beta_SM1550 = -19.118488
beta_DCF = 73.924818
beta_DSF = 5.098263
dispersion1 = np.exp(1j*beta_SM1550/2*(omiga_all**2)*L1)
dispersion2 = np.exp(1j*beta_DCF/2*(omiga_all**2)*L2)
dispersion3 = np.exp(1j*beta_DSF/2*(omiga_all**2)*L3)
aa = np.random.uniform(0, 1, N)
bb = np.random.uniform(0, 1, N)
u = 0.1*(aa+1j*bb)
Ep = np.trapz(t, np.abs(u)**2)
Puu = Ep/T_AM
u_begin = u
# pulse evolution
dz = 0.0001
tol = 1e-6
N_trip = 500
spec_z = np.empty((N_trip, N))
u_z = np.empty((N_trip, N))
for _ in trange(N_trip):
    u1 = u
    u = u*T_intense*T_phase*0.5
    u_temp1 = fft.fftshift(fft.ifft(fft.fftshift(u)))*dispersion1
    u = fft.fftshift(fft.fft(fft.fftshift(u_temp1)))

    u, nf, Plotdata_EDF = IP_CQEM_FD(u,dt,dz,EDF_1550,fo,tol,1,1) # across the EDF
   
    u_temp2 = dispersion2*fft.fftshift(fft.ifft(fft.fftshift(u)))
    u = fft.fftshift(fft.fft(fft.fftshift(u_temp2)))

    u_temp3 = dispersion3*fft.fftshift(fft.ifft(fft.fftshift(u)))
    u = fft.fftshift(fft.fft(fft.fftshift(u_temp3)))

    uout, ur = opti_coupler(u, 0, opti_oc.rho_out)
    u = ur
    u_OC = u
    u = filter_gauss(u, opti_filter.f3dB, opti_filter.fc, opti_filter.n, fo, df)
    
    spec = abs(Et2Ef(uout, T)**2)
    specnorm = spec/lambda_1**2
    spec_z[_:] = specnorm
    u_z[_:] = np.real(uout)
save_data = pd.DataFrame({'t': t, 'u': u, 'lambda': lambda_all, 'u_begin': u_begin})
save_data.to_csv('d:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/data/pulse_evolution_result.csv')
plt.figure(1)
# plt.plot(t, np.abs(u_begin), 'r', label='before')
plt.plot(t, np.abs(u), 'b', label='after')
plt.legend()
plt.xlabel('t(ps)')
plt.ylabel('P')
plt.title('Pulse evolution result')
plt.grid()
plt.xlim(-30, 30)
plt.savefig('d:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/data/pulse_evolution_result.png')


plt.figure(2)
plt.plot(lambda_all, 10*np.log10(
    np.abs(fft.fftshift(fft.fft(u)))**2*1e3))
plt.xlabel('lambda(nm)')
plt.ylabel('P(dBm)')
plt.title('Spectrum evolution result')
plt.savefig('d:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/data/Spectrum_evolution_result.png')
plt.grid()
# plt.xlim(1540, 1560)
# plt.ylim(-250, 0)

fig3, ax3 = plt.subplots(1, 2)
fig3.set_size_inches(10, 5.5)
c1 = ax3[0].pcolor(t, np.arange(N_trip), np.abs(u_z)**2, cmap='jet')
ax3[0].set_title('The Output Signal Field Evolution')
ax3[0].set_xlabel('Time/ps')
ax3[0].set_ylabel('Trips')
fig3.colorbar(c1, ax=ax3[0])
c1 = ax3[1].pcolor(c/(f+fo), np.arange(N_trip), spec_z, cmap='jet')
ax3[1].set_title('The Output Spectrum Field Evolution')
ax3[1].set_xlabel('Wavelength/nm')
ax3[1].set_ylabel('Trips')
fig3.colorbar(c1, ax=ax3[1])

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

plt.show()


