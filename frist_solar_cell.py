from solcore import material
from solcore.structure import Layer, Junction, TunnelJunction
from solcore.solar_cell import SolarCell
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource
import solcore.poisson_drift_diffusion as PDD

import numpy as np
import matplotlib.pyplot as plt

T = 298
wl = np.linspace(350, 1050, 301) * 1e-9


n_GaAs = material("GaAs")(T=T, Nd=1e24, electron_diffusion_length= 2e-6,hole_diffusion_length=0.1e-6,
                          hole_mobility=5e-2,
                          electron_mobility = 3.4e-3,
                            )
p_GaAs = material("GaAs")(T=T, Na=8e22, hole_diffusion_length= 0.1e-6 , electron_diffusion_length= 2e-6,
                          hole_mobility=3.4e-3,
                          electron_mobility=5e-2,
                          )
GaAs_junction = Junction(
    [
        Layer(width=150e-9, material=p_GaAs, role="Emitter"),
        Layer(width=2000e-9, material=n_GaAs, role="Base"),
    ],
    T=T,
    kind="PDD",
)
my_solar_cell = SolarCell(
    [
        GaAs_junction,
    ],
    T=T,
    # substrate=p_GaAs,
)

light_source = LightSource(
    source_type="standard",
    version="AM1.5g",
    x=wl,
    output_units="photon_flux_per_m",
    concentration=1,
)

solar_cell_solver(my_solar_cell,"qe",user_options={
        "light_source": light_source,
        "wavelength": wl,
        "optics_method": "TMM",
    },)

vint = np.linspace(-3.5, 4, 600)
V = np.linspace(0, 3.5, 300)


fig3, axIV = plt.subplots(1, 1, figsize=(6, 4))
solar_cell_solver(
    my_solar_cell,#วัสดุ
    "iv",#วัดIV
    user_options={
        "light_source": light_source,
        "wavelength": wl,
        "optics_method": None,
        "light_iv": True,#ให่IV ออกมา
        "mpp": True,#ให้Mppออกมา
        "voltages": V,#ให้Vออกมา
        "internal_voltages": vint,
    },
)
axIV.plot(V, my_solar_cell.iv["IV"][1] /my_solar_cell.iv["Isc"], )#เป็นการวัดลักษณะIVเมิอแสงเปลี่ยนไป

axIV.legend(loc="lower left", frameon=False)
axIV.set_ylim(0, 1.1)
axIV.set_xlim(0, 3.5)
axIV.set_xlabel("Voltage (V)")
axIV.set_ylabel("Normalized current (-)")

plt.tight_layout()

fig, ax1 = plt.subplots(1, 1, figsize=(6, 4))

ax1.plot(wl * 1e9, my_solar_cell.absorbed, "k", label="Total Absorbed")
ax1.legend(loc="upper right", frameon=False)
ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel("EQE")
ax1.set_ylim(0, 1.1)
ax1.set_xlim(350, 1150)

plt.tight_layout()

plt.show()
