import sys

import numpy as np
import scipy as sc

sys.path.append('D:\Onedrive\OneDrive - CQU\Code\Python\Activate_Mode_Lock\methods')
from methods_simplify.Optical_functions import *


def IP_CQEM_FD(u0, dt, dz, mod, fo, tol, dplot, quiet):
    # sourcery skip: low-code-quality, merge-else-if-into-elif
    nt = len(u0)
    w = np.fft.fftshift(2*np.pi*np.linspace(-nt/2, nt/2-1, nt)/(dt*nt)) # angular frequency
    t = np.linspace(0, (nt-1)*dt, nt)                   # time
    t = t - t[nt//2]                                    # time centered at zero
    u1 = u0
    hrw,fr = Raman_response_w(t,mod)                    # Raman response

    # preparation before the IPM
    ufft = np.fft.fft(u0)                               # FFT of the input signal
    # ufft = np.fft.fftshift(ufft)                        # shift the FFT to the center
    propagedlength = 0                                  # the length of the fiber that has been propaged
    nf = 1                                              

    if 'gssdB' in mod.__dict__.keys():
        gain_w = filter_lorenz_tf(u1,mod.fbw,mod.fc,fo,1/dt/nt)
        alpha_0 = mod.alpha

    if dplot == 1:
        z_all = np.array([])
        ufft_z = np.array([])
        u_z = np.array([])

    if quiet == 0:
        print('IPM:')
    alpha_new = mod.alpha
    while propagedlength < mod.length:
        if (dz+propagedlength) > mod.length:
            dz = mod.length - propagedlength
        if 'gssdB' in mod.__dict__.keys():
            Pin0 = (sum(u1*np.conj(u1))/nt).real
            gain = gain_saturated(Pin0, mod.gssdB, mod.PsatdBm)*gain_w
            alpha_new = alpha_0 - gain
        if type(alpha_new) == np.ndarray:
            LOP = Linearoperator_w(alpha_new, mod.betaw, w)
        else:
            LOP = Linearoperator_w_c(alpha_new, mod.betaw, w)


        PhotonN = sum(abs(ufft)**2/(w+2*np.pi*fo))
        if type(alpha_new) == np.ndarray:
            Photon_z = sum(np.exp(-dz*np.fft.fftshift(alpha_new))*abs(ufft)**2/(w+2*np.pi*fo))
        else:
            Photon_z = sum(np.exp(-dz*alpha_new)*abs(ufft)**2/(w+2*np.pi*fo))

        halfstep = np.exp(LOP*dz/2)
        uip = halfstep*ufft
        k1 = halfstep*dz*NonLinearoperator_w(u1, mod.gamma,w, fo, fr, hrw, dt,mod)

        uhalf2 = np.fft.ifft(uip + k1/2)
        k2 = dz*NonLinearoperator_w(uhalf2, mod.gamma, w, fo, fr, hrw, dt, mod)

        uhalf3 = np.fft.ifft(uip + k2/2)
        k3 = dz*NonLinearoperator_w(uhalf3, mod.gamma, w, fo, fr, hrw, dt, mod)

        uhalf4 = np.fft.ifft(halfstep*(uip + k3))
        k4 = dz*NonLinearoperator_w(uhalf4, mod.gamma, w, fo, fr, hrw, dt, mod)

        uaux = halfstep*(uip + k1/6 + k2/3 + k3/3) + k4/6;    

        propagedlength += dz

        if quiet == 0:
            print('propagedlength = ', propagedlength*100/mod.length, 'm')
        error1 = abs(sum(abs(uaux)**2/(w+2*np.pi*fo)) - Photon_z)/Photon_z
        if error1 > 2*tol:
            propagedlength -= dz
            dz = dz/2
        else:
            ufft = uaux
            u1 = np.fft.ifft(ufft)
            if error1 > tol:
                dz = dz/(2**0.2)
            else:
                if error1 < tol/2:
                    dz = dz*2**0.2
            if dplot == 1:
                z_all = np.append(z_all, propagedlength)
                ufft_z = np.append(ufft_z, abs(np.fft.fftshift(ufft)))
                u_z = np.append(u_z, u1)
        nf += 16
    class data():
        def __init__(self):
            self.z = 0
            self.ufft = 0
            self.u = 0
    Plotdata = data()
    if dplot == 1:
        Plotdata.z = z_all
        Plotdata.ufft = ufft_z
        Plotdata.u = abs(u_z)
    else:
        Plotdata = 0
    # print(Plotdata.u.shape, Plotdata.ufft.shape)
    # Plotdata.u = Plotdata.u.reshape((4096, 7))
    # Plotdata.ufft = np.reshape(Plotdata.ufft, (4096, 7))
    return u1, nf, Plotdata