# from solcore.absorption_calculator import calculate_absorption_profile
from solcore import material
from solcore.structure import Layer, Junction, TunnelJunction
from solcore.solar_cell import SolarCell
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource
# import solcore.poisson_drift_diffusion as PDD
import numpy as np
import matplotlib.pyplot as plt
import pickle
MgF2 = material("MgF2")() #1.3917
ZnS = material("ZnScub")() #2.4199
Y2O3 = material('Y2O3')() #1.9822
Ta2O5 = material('TAOX1')() #2.09
n_GaAs = material("GaAs")(T=300, Nd=1e24) #3.66
p_GaAs = material("GaAs")(T=300, Na=8e22)
wl = np.linspace(300, 1000, 200) *1e-9#wavelength
light_source = LightSource(
    source_type="standard",
    version="AM1.5g",
    x=wl,
    output_units="photon_flux_per_m",
    concentration=1,
)

T = 300
window_bottom = material("GaInP")(T=T, Nd=5e24, In=0.49)
# n_GaAs = material("GaAs")(T=300, Nd=1e18) #3.66
# p_GaAs = material("GaAs")(T=300, Na=8e16)
n_GaAs = material("GaAs")(T=T, Nd=1e24)
p_GaAs = material("GaAs")(T=T, Na=8e22)

bsf_bottom = material("GaInP")(T=T, Na=5e24, In=0.49)

vint = np.linspace(-3.5, 4, 600)
V = np.linspace(-3.5, 0, 300)

GaAs_junction = Junction([
    # Layer(width=10e-9, material=window_bottom, role="Window"),  #In width Nd
    Layer(width=150e-9, material=n_GaAs, role="Emitter"),       #Nd width
    Layer(width=2000e-9, material=p_GaAs, role="Base"),         #Na width
    # Layer(width=200e-9, material=bsf_bottom, role="BSF"),       #Na In
],T=300,kind="PDD",sn=1e6,sp=1e6)


my_solar_cell = SolarCell([
    Layer(width= 90e-9 , material=MgF2, role = 'ARC1'),
    Layer(width= 50e-9 , material=ZnS, role = 'ARC2'),
    GaAs_junction]
    ,T=300
    ,substrate=p_GaAs,)

solar_cell_solver(my_solar_cell,"qe",
    user_options={"light_source": light_source,
        "wavelength": wl,
        "optics_method": "TMM",
    },)
solar_cell_solver(my_solar_cell,"iv",
               user_options={
                   "light_source": light_source,
                   "wavelength": wl,
                   "optics_method": None,
                   "light_iv": True,
                   "mpp": True,
                   "voltages": V,
                   "internal_voltages": vint,},)
print(my_solar_cell.iv['Pmpp'])
print(my_solar_cell.iv['Isc'])
print(my_solar_cell.iv['Voc'])
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.25, 4))
for j in my_solar_cell.junction_indices:#junctionคือหยั่ง
    zz = (
        my_solar_cell[j].short_circuit_data.Bandstructure["x"] + my_solar_cell[j].offset)
    n = my_solar_cell[j].short_circuit_data.Bandstructure["n"]
    p = my_solar_cell[j].short_circuit_data.Bandstructure["p"]
    ax1.semilogy(zz * 1e9, n, "b")#อันนี้น่าจะเป็นการระบุcarrier densityของตอนฉายแสงของe และhole
    ax1.semilogy(zz * 1e9, p, "r")

    zz = my_solar_cell[j].equilibrium_data.Bandstructure["x"] + my_solar_cell[j].offset
    n = my_solar_cell[j].equilibrium_data.Bandstructure["n"]
    p = my_solar_cell[j].equilibrium_data.Bandstructure["p"]
    ax1.semilogy(zz * 1e9, n, "b--")#อันนี้น่าจะเป็นการระบุcarrier densityของตอนปิดแสงของe และhole
    ax1.semilogy(zz * 1e9, p, "r--")

ax1.set_xlabel("Position (nm)")
ax1.set_ylabel("Carrier density (m$^{-3}$)")
plt.tight_layout()
labels = ["EQE GaAs"]
colours = ["b", "r"]
for i, j in enumerate(my_solar_cell.junction_indices):#เป็นการพล็อตของeqeของแต่ละอันของGaAs และGaInP
    ax2.plot(wl * 1e9, my_solar_cell[j].eqe(wl), colours[i], label=labels[i])

ax2.plot(wl * 1e9, my_solar_cell.absorbed, "k", label="Total Absorbed")

ax2.legend(loc="upper right", frameon=False)
ax2.set_xlabel("Wavelength (nm)")
ax2.set_ylabel("EQE")
ax2.set_ylim(0, 1.1)
ax2.set_xlim(350, 1150)
plt.tight_layout()
plt.show()
with open('data_solar_cell.pickle', 'wb') as out:
    pickle.dump(my_solar_cell, out)
# print(my_solar_cell.absorbed)
