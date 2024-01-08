import numpy as np
import matplotlib.pyplot as plt

from solcore.light_source import LightSource
import solcore.light_source.smarts
from solcore import config

config.spice('C:/Users/usuwa/Desktop/seneir project/simulation/solcore5-develop/solcore5-develop/solcore/spice')
# config.smarts('C:/Users/usuwa/PycharmProjects/smarts-295-pc/SMARTS_295_PC')
print(config)
# The wavelength range of the spectra
wl = np.linspace(300, 3000, 200)

# Now different types of light sources can be defined
gauss = LightSource(source_type='laser', x=wl, center=800, linewidth=50, power=200)
bb = LightSource(source_type='black body', x=wl, T=5800, entendue='Sun')
am15g = LightSource(source_type='standard', x=wl, version='AM1.5g')
# smarts = LightSource(source_type='SMARTS', x=wl)
spectral = LightSource(source_type='SPECTRAL2', x=wl)

# Plot comparing the different light sources
plt.figure(1)
plt.plot(*gauss.spectrum(), label='Gauss')
plt.plot(*bb.spectrum(), label='Black body')
plt.plot(*am15g.spectrum(), label='AM1.5G')
# plt.plot(*smarts.spectrum(), label='SMARTS')
plt.plot(*spectral.spectrum(), label='SPECTRAL2')

plt.xlim(300, 3000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Power density (Wm$^{-2}$nm$^{-1}$)')
plt.tight_layout()
plt.legend()

# Plot comparing the spectra calculated with SMARTS at different hours of the day
plt.figure(2)
# for h in range(8, 20):
#     plt.plot(*smarts.spectrum(HOUR=h), label='{} h'.format(h))

plt.xlim(300, 3000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Power density (Wm$^{-2}$nm$^{-1}$)')
plt.tight_layout()
plt.legend()
plt.show()
