import math


class Modulator():
    def __init__(self, name, Vpi, Vbias, phase) -> None:
        # the input parameters stand for the types, the Vpi, the bias voltage, the phase
        self.name = name
        self.Vpi = Vpi
        self.Vbias = Vbias
        self.phase = phase
        self.Va = self.Vpi/2
        self.Gamb = math.pi*self.Vbias/self.Vpi # the gain of the modulator
        self.Gamm = math.pi*self.Va/self.Vpi