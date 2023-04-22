import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt

fs = 1e9
Nt = 2**10
dt = 1/fs
t = np.linspace(0, Nt*dt, Nt)
df = fs/Nt
f = np.linspace(-fs/2, fs/2-df, Nt)
fsignal = 4e8
signal = np.sin(2*np.pi*fsignal*t)
freq = np.fft.fftshift(np.fft.fft(signal))
plt.plot(f, np.abs(freq))
plt.show()
