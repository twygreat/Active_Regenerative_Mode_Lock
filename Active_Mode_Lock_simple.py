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
freq_all = f + fo
lambda_all = c/freq_all
omiga_all = 2*np.pi/T*np.linspace(-N/2, N/2-1, N)

# AM parameters
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

# fiber parameters

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
EDF_1550.betaw = [0, 0, 1.7, 0.171]

g0 = 2
Psat = 50/1000
L1 = 10
L2 = 0.5
L3 = 20
beta_SM1550 = -1.01
beta_DCF = 4.49
beta_DSF = 0.34
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
for _ in trange(1000):
    u1 = u
    u = u*T_intense*T_phase*0.5
    u_temp1 = fft.fftshift(fft.ifft(fft.fftshift(u)))*dispersion1
    u = fft.fftshift(fft.fft(fft.fftshift(u_temp1)))

    u, nf, Plotdata_EDF = IP_CQEM_FD(u,dt,dz,EDF_1550,fo,tol,1,1) # across the EDF
   
    # u_temp2 = dispersion2*fft.fftshift(fft.ifft(fft.fftshift(u)))
    # u = fft.fftshift(fft.fft(fft.fftshift(u_temp2)))

    u_temp3 = dispersion3*fft.fftshift(fft.ifft(fft.fftshift(u)))
    u = fft.fftshift(fft.fft(fft.fftshift(u_temp3)))

save_data = pd.DataFrame({'t': t, 'u': u, 'lambda': lambda_all, 'u_begin': u_begin})
save_data.to_csv('d:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/z_data/pulse_evolution_result.csv')
plt.figure(1)
# plt.plot(t, np.abs(u_begin), 'r', label='before')
plt.plot(t, 10*np.log10(np.abs(u)), 'b', label='after')
plt.legend()
plt.xlabel('t(ps)')
plt.ylabel('P')
plt.title('Pulse evolution result')
plt.grid()
plt.savefig('d:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/z_data/pulse_evolution_result.png')


plt.figure(2)
plt.plot(lambda_all, 10*np.log10(
    np.abs(fft.fftshift(fft.ifft(u)))**2*1e3))
plt.xlabel('lambda(nm)')
plt.ylabel('P(dBm)')
plt.title('Spectrum evolution result')
plt.savefig('d:/Onedrive/OneDrive - CQU/Code/Python/Active_Regenerative_Mode_Lock/z_data/Spectrum_evolution_result.png')
plt.grid()
# plt.xlim(1540, 1560)


plt.show()


