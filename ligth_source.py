import numpy as np
import matplotlib.pyplot as plt

from solcore.light_source import LightSource

#----------------------------
#def
h = 4.14*10**(-15)
c = 3*10**8
q = 1.6*10**-19

def ligthwave2ev(wave):
    wave = np.divide(wave,10**9)
    output = np.divide(h*c,wave)
    return output

wl = np.linspace(300,3000,1000)
AM15 = LightSource(source_type="standard", x=wl, version="AM1.5g", concentration=1,)
AM0  = LightSource(source_type="standard", x=wl, version="AM0", concentration=1,)
BB   = LightSource(source_type="black body",x=wl,  T=5800,entendue='Sun')
print(ligthwave2ev(wl))
plt.figure(1)
plt.plot(*AM15.spectrum(), label='AM 1.5')
# plt.plot(*AM0(), label='AM 0')
plt.plot(*BB.spectrum(), label='Black Body,Ts=5800K')
plt.xlim(300,3000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Power density (Wm$^{-2}$nm$^{-1}$)')
plt.tight_layout()
plt.legend()
plt.show()