import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, trapz
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

wl   = np.linspace(300,3000,1000)
AM15 = LightSource(source_type="standard",  x=wl, version="AM1.5g", concentration=1,)
AM0  = LightSource(source_type="standard",  x=wl, version="AM0", concentration=1,)
BB   = LightSource(source_type="black body",x=wl, T=5800, entendue='Sun')
BB_ev= LightSource(source_type='black body',x=wl, T=5800, entendue='Sun', output_units="power_density_per_ev")
# print(ligthwave2ev(wl))
fig, ax_nm = plt.subplots(1,1,figsize=(8, 6))
light_source_measure = LightSource(
    source_type="standard",
    version="AM1.5g",
    output_units='power_density_per_nm',
    x=wl,
    concentration=1,
)
spectrum = light_source_measure.spectrum()
# Calculate the power density in 1 cm^2
area_cm2 = 1  # cm^2
area_m2 = area_cm2 * (1e-2)**2  # Convert cm^2 to m^2

print(spectrum)
# Integrate to get the power in watts
power_watts = trapz(spectrum, wl) * area_m2
print("Power in 1 cm^2 for AM1.5G spectrum:", power_watts, "W")

spectrum = AM15.spectrum()[1]

# Calculate the power density in 1 cm^2
area_cm2 = 1  # cm^2
area_m2 = area_cm2 * (1e-2)**2  # Convert cm^2 to m^2

# Integrate to get the power in watts
power_watts = trapz(spectrum, wl) * area_m2

print("Power in 1 cm^2 for AM1.5G spectrum:", power_watts*1000, "mW")
ax_nm.plot(*AM15.spectrum(), label='AM 1.5')
# ax_nm.plot(*AM0(), label='AM 0')
ax_nm.plot(*BB.spectrum(), label='Black Body,Ts=5800K')
ax_nm.set_xlim(300,3000)
ax_nm.set_xlabel('Wavelength (nm)')
ax_nm.set_ylabel('Power density (Wm$^{-2}$nm$^{-1}$)')
ax_nm.legend()

# print(*BB_ev.spectrum())
# fig2, axev = plt.subplot(1, 1, figsize=(6, 4))
# axev.plot(*BB_ev.spectrum(), label="Black Body,Ts=5800K")
plt.tight_layout()
plt.show()