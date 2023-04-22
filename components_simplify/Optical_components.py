"""All optical components are defined here"""
class fiber():
    n = 1.45 # the refractive index of the fiber
    c = 3e8 # the speed of light
    lambda_pulse = 1550

    def __init__(self, id, length, diameter, aeff, gamma, alpha, beta2, n2):
        self.id = id # the id of the fiber
        self.length = length # the length of the fiber 
        self.diameter = diameter # the diameter of the fiber
        self.aeff = aeff # the effective area of the fiber
        self.gamma = gamma # the nonlinear coefficient of the fiber
        self.alpha = alpha # the attenuation coefficient of the fiber
        self.beta2 = beta2 # the dispersion coefficient of the fiber
        self.n2 = n2 # Kerr coefficient (10^-20*m^2/W)

    def c_beta2(self, l):
        return l*self.beta2
    

class NarrowFilter():
    def __init__(self, Gwidh) -> None:
        self.Gwidh = Gwidh


class OpticalFilter():
    c = 299792.458
    def __init__(self, lamda_c, lamda_bw, n):
        self.lamda_c = lamda_c
        self.lamda_bw = lamda_bw
        self.fc = self.c/self.lamda_c
        self.f3dB = self.c/(self.lamda_c)**2*self.lamda_bw
        self.n = n


class OptiOc():
    def __init__(self, rho_out) -> None:
        self.rho_out = rho_out
        self.out = 1 - self.rho_out


class Absorber():
    def __init__(self, name,q_0, q_1, P_0,Peak_Power):
        self.name = name
        self.q_0 = q_0
        self.q_1 = q_1
        self.P_0 = P_0
        self.Peak_Power = Peak_Power