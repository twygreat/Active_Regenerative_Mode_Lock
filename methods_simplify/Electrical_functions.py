import numpy as np


def Modulation_way_1(f_Thz, Gamb, Gamm, t, Phase_I):
    return 1+np.exp(1j*(Gamb+Gamm*np.cos(2*np.pi*f_Thz*t+Phase_I)))

def Modulation_way_2(f_Thz, Phase_P, t, Vpi_P, Vp):
    Phase_out = np.pi*Vp*np.cos(2*np.pi*f_Thz*t+Phase_P)/Vpi_P
    return np.exp(1j*Phase_out)

def Modulation_way_simple(u, Vbias, Vi, Vp, N, f, Vin, t):
    b = Vbias/Vi
    m_i = Vin/Vi
    m_p = Vin/Vp
    w_i = 2*np.pi*f/1000
    nn = min(2*round(500/f*(N/1000)), N)
    piecewise_1 = np.zeros(int((N-nn)/2))
    piecewise_2 = np.ones(int(nn))
    piecewise_3 = np.zeros(int((N-nn)/2))
    piecewise = np.concatenate((piecewise_1, piecewise_2, piecewise_3))
    T_intense = 0.5*(1+np.sin(np.pi/2*(b+m_i*np.cos(w_i*t))))*piecewise
    T_phase = np.exp(1j*(np.pi*m_p*np.cos(w_i*t)))*piecewise
    u = u*T_intense*T_phase*0.5
    return u
