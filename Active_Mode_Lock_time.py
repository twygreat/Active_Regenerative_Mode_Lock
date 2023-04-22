import numpy as np
import matplotlib.pyplot as plt
import scipy.fftpack as fft
from tqdm import trange

def shah(w, dw, n):
    dw_seq = np.zeros(n)
    m = round(dw/w*n)
    for i in range(0, n, m):
        dw_seq[i] = 1
    return dw_seq

def fguassian(FWHM, N, dt):
    t = np.linspace(-N/2*dt, (N/2-1)*dt, N)
    T = FWHM /2/np.sqrt(np.log(2)/2)
    Eenv = np.exp(-(t/T)**2)
    return Eenv, t

def Efft(Et, t):
    nt = len(Et)
    T = t[-1] - t[0]
    dft = 1/T
    n = np.arange(-nt/2, nt/2)
    ft = n*dft
    Eft = fft.fftshift(fft.fft(Et))
    return Eft, ft

def Eifft(Eft, ft):
    nt = len(Eft)
    Eft = fft.ifftshift(Eft)
    Et = fft.ifft(Eft, nt)
    return Et


def main():
    nt = 2**14
    c = 3e8
    wavelength = 1550e-9
    l_all = 0.025 # the length of the resonator
    d_mf = c/1.45/l_all
    rpt = 1/d_mf*1e9
    print(f'基模频率={d_mf/1e9} GHz')
    print(f'周期={rpt}ns')
    t_FWHM = 220e-15
    dt = 50e-15
    d_WB = 0.441/t_FWHM*wavelength**2/c
    Et, t = fguassian(t_FWHM, nt, dt)
    It = Et**2
    It = It/max(It)
    plt.figure()
    plt.plot(t*1e15, It)
    plt.xlabel('t (fs)')
    plt.ylabel('I (a.u.)')
    plt.title('time domain')
    plt.xlim(-2000, 2000)
    plt.grid()

    Eft, ft = Efft(Et, t)
    ang_w = 2*np.pi*(c/wavelength+ft)
    wlgth = c/ang_w*2*np.pi
    I_w = Eft*np.conj(Eft)*c/wlgth**2
    I_w = I_w/max(I_w)
    plt.figure()
    plt.plot(wlgth*1e9, abs(I_w))
    plt.xlabel('wavelength (nm)')
    plt.ylabel('I (a.u.)')
    plt.title('frequency domain')
    plt.grid()

    delta_F = ft[-1]-ft[0]
    delta_Seq = shah(delta_F, d_mf, nt)
    Eft = Eft*delta_Seq
    I_wlgth = Eft*np.conj(Eft)*c/wlgth**2
    I_wlgth = I_wlgth/max(I_wlgth)
    plt.figure()
    plt.plot(wlgth*1e9, abs(I_wlgth))
    plt.xlabel('wavelength (nm)')
    plt.ylabel('I (a.u.)')
    plt.title('frequency domain actual')
    plt.grid()

    Et = Eifft(Eft, ft)
    It = Et*np.conj(Et)
    plt.figure()
    plt.plot(t*1e9, abs(It))
    plt.xlabel('t (ns)')
    plt.ylabel('I (a.u.)')
    plt.title('time domain actual')
    plt.grid()

    ns_phase = np.random.uniform(0, 2*np.pi, nt)
    Eft = Eft*np.exp(1j*ns_phase)
    Et = Eifft(Eft, ft)
    It = Et*np.conj(Et)
    It = It/max(It)
    plt.figure()
    plt.plot(t*1e9, abs(It))
    plt.xlabel('t (ns)')
    plt.ylabel('I (a.u.)')
    plt.ylim(-1, 6)
    plt.title('time domain origin')


    plt.show()


if __name__ == '__main__':
    main()
