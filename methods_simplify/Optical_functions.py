"""All optical functions are defined here"""
import math

import numpy as np
import scipy as sp


def Et2Ef(Et,twidow):
    Ef = sp.fft.fftshift(sp.fft.ifft(sp.fft.fftshift(Et)))*twidow
    return Ef

def filter_gauss(ui, f3dB, fc, n, fo, df):
    Ui = np.fft.fft(ui)
    N = Ui.shape[0]
    f = np.linspace(-N/2, N/2-1, N) * df + fo
    f = np.fft.fftshift(f)
    w = 2*np.pi*f
    Tf = np.exp(-np.log(math.sqrt(2))*(2/f3dB*(f-fc))**(2*n))
    return np.fft.ifft(Ui*Tf)

def filter_lorenz_tf(ui,fbw,fc,fo,df):
    N = ui.shape[0]
    f = np.linspace(-N/2,N/2-1,N)*df + fo
    w = 2*np.pi*f
    tf = fbw/2/np.pi/((f-fc)**2+(fbw/2)**2)
    tf = tf/np.max(tf)
    return tf

def gain_saturated(Pin, gssdB, PsatdBm):
    """
    This function is used to calculate the gain of the fiber
    """
    gss = 10**(gssdB/10)
    Psat = 10**(PsatdBm/10)/1000
    return gss/(1+Pin/Psat)

def Linearoperator_w_c(alpha, betaw, w):
    LOP= -alpha/2
    if len(betaw) == len(w):
        LOP = LOP - 1j*betaw
        LOP = np.fft.fftshift(LOP)
    else:
        for i in range(len(betaw)):
            LOP = LOP - 1j*betaw[i]*w**i/np.math.factorial(i)
    return LOP

def Linearoperator_w(alpha, betaw, w):
    LOP= -np.fft.fftshift(alpha/2)
    if len(betaw) == len(w):
        LOP = LOP - 1j*betaw
        LOP = np.fft.fftshift(LOP)
    else:
        for i in range(len(betaw)):
            LOP = LOP - 1j*betaw[i]*w**i/np.math.factorial(i)
    return LOP

def NonLinearoperator_w(u_t, gamma, w, fo, fr, hrw, dt, mod):
    if 'ssp' not in mod.__dict__.keys() or mod.ssp == 1:
        NLOP = -1j*gamma*(1+w/(2*np.pi*fo))*np.fft.fft((1-fr)*u_t*abs(u_t)**2)\
            +fr*dt*u_t*np.fft.ifft(hrw*np.fft.fft((u_t)**2))
    else:
        NLOP = -1j*gamma*np.fft.fft((1-fr)*u_t*abs(u_t)**2)\
            +fr*dt*u_t*np.fft.ifft(hrw*np.fft.fft((u_t)**2))
        
    return NLOP

def nonlinerphaseshift(input_arg1, input_arg2):
    # sourcery skip: inline-immediately-returned-variable
    data_1 = input_arg1
    Peak_Power = max(abs(data_1.u)**2)
    dz = data_1.z[1]-data_1.z[0]
    out_arg = sum(input_arg2.gamma*Peak_Power*dz)
    return out_arg

def opti_coupler(u1i, u2i, rho):
    if rho > 1:
        rho = 1
    elif rho < 0:
        rho = 0
    u1o = np.sqrt(rho)*u1i + 1j*np.sqrt(1-rho)*u2i
    u2o = np.sqrt(rho)*u2i + 1j*np.sqrt(1-rho)*u1i
    return u1o, u2o

def Raman_response_w(t, mod):
    if 'raman' in mod.__dict__.keys() and mod.raman == 0:
        hrw = 0
        fr = 0
    else:
        # raman parameters
        t1 = 12.2e-3
        t2 = 32e-3
        tb = 96e-3
        fc = 0.04
        fb = 0.21
        fa = 1-fc-fb
        fr = 0.245
        tres = t-t[0]
        ha = ((t1**2+t2**2)/(t1*t2**2))*np.exp(-tres/t2)*np.sin(tres/t1)
        hb = ((2*tb-tres)/tb**2)*np.exp(-tres/tb)
        hr = (fa+fc)*ha+fb*hb
        hrw = np.fft.fftshift(np.fft.fft(hr))
    return hrw, fr

def Saturation(At, q_0, q_1, P_0):
    return np.sqrt(1-(q_0/(1+At*np.conj(At)/P_0))-q_1)

def Aoppm(T, w):  # sourcery skip: inline-immediately-returned-variable
    h = 6.626e-34
    N = len(T)
    dT = T[1]-T[0]
    Twindow = T[-1]-T[0]
    w0 = 2*np.pi*np.linspace(-N/2,N/2-1,N)/(N*dT)
    nu = (w0+w)/(2*np.pi)
    phi = 2*np.pi*np.random.rand(N)
    Aoppm_f = np.sqrt(h*nu*Twindow)*np.exp(1j*phi)
    answer = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(Aoppm_f)))
    return answer

def fwhm(x):
    peak, index_peak = max(x), np.argmax(x)
    half_peak = peak/2
    x_1 = x[:index_peak]
    x_2 = x[index_peak:]
    I_1 = np.argwhere(x_1 <= half_peak)
    I_2 = np.argwhere(x_2 <= half_peak)
    I_1 = int(I_1[-1])
    I_2 = int(I_2[0])
    width = index_peak-I_1+I_2
    I_l = I_1
    I_r = index_peak + I_2
    return width, I_l, I_r