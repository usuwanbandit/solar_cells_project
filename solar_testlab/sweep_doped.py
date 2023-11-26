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
window_bottom = material("GaInP")(T=T, Nd=5e24, In=0.49)
n_GaAs = material("GaAs")(T=T, Nd=1e24)
p_GaAs = material("GaAs")(T=T, Na=8e22)
bsf_bottom = material("GaInP")(T=T, Na=5e24, In=0.49)
MgF2 = material("MgF2")() #1.3917
ZnS = material("ZnScub")() #2.4199
light_source = LightSource(
    source_type="standard",
    version="AM1.5g",
    x=wl,
    output_units="photon_flux_per_m",
    concentration=1,
)
vint = np.linspace(-2, 1, 300)
V = np.linspace(-2, 0, 150)
def solver_GaAs_for_sweep(Nd, Na):
    n_GaAs = material("GaAs")(T=T, Nd=Nd)
    p_GaAs = material("GaAs")(T=T, Na=Na)
    GaAs_junction = Junction(
        [
            Layer(width=10e-9, material=window_bottom, role="Window"),  # In width Nd
            Layer(width=150e-9, material=n_GaAs, role="Emitter"),  # Nd width
            Layer(width=2000e-9, material=p_GaAs, role="Base"),  # Na width
            Layer(width=200e-9, material=bsf_bottom, role="BSF"),  # Na In
        ],
        T=T,
        kind="PDD")
    my_solar_cell = SolarCell([
        Layer(width= 90e-9 , material=MgF2, role = 'ARC1'),
        Layer(width= 50e-9 , material=ZnS, role = 'ARC2'),
        GaAs_junction
    ], T=T, substrate=p_GaAs, )
    solar_cell_solver(my_solar_cell, "iv",
                      user_options={
                          "light_source": light_source,
                          "wavelength": wl,
                          "optics_method": 'TMM',
                          "light_iv": True,
                          "mpp": True,
                          "voltages": V,
                          "internal_voltages": vint,
                      }, )
    return my_solar_cell

# print(solver_GaAs_for_sweep(1e28,2e22).iv["Pmpp"])
def sweep():
    doped_emitter_num = 14
    doped_base_num = 14

    doped_emitter_con = np.power(10, np.linspace(14, 28, doped_emitter_num))
    doped_base_con = np.power(10, np.linspace(14, 28, doped_base_num))
    # vint = np.linspace(-2, 1, 300)
    # V = np.linspace(-2, 0, 150)
    size = (doped_emitter_num, doped_base_num)
    isc_np = np.zeros(size)
    voc_np = np.zeros(size)
    FF_np = np.zeros(size)
    pmpp_np = np.zeros(size)
    # allI_np = np.zero((width_con, doped_con))
    index = 0
    N = doped_emitter_num * doped_base_num
    # input('press any buttom to start')
    for i, doped_emitter in enumerate(doped_emitter_con):
        for j, doped_base in enumerate(doped_base_con):
            my_solar_cell = solver_GaAs_for_sweep(doped_emitter,doped_base)
            isc_np[i,j] =  my_solar_cell.iv["Isc"]
            voc_np[i,j] = my_solar_cell.iv["Voc"]
            FF_np[i, j]  = my_solar_cell.iv["FF"]
            pmpp_np[i, j] = my_solar_cell.iv["Pmpp"]
            # allI_np[i, j] =my_solar_cell.iv["IV"][1]
            # emitter_width(ii)
            index += 1
            print('==============================\n')
            print(int(index / N * 100), "%\n")
            print(f'complete {index} out of {N}')
            print('==============================\n')
    with open('../data.npy', 'wb') as fout:
        np.save(fout, isc_np)
        np.save(fout, voc_np)
        np.save(fout, FF_np)
        np.save(fout, pmpp_np)
from numpy import ma
from matplotlib import cm, ticker
with open('../data.npy', 'rb') as fin:
    isc_np = np.load(fin)
    voc_np= np.load(fin)
    FF_np= np.load(fin)
    pmpp_np = np.load(fin)
doped_emitter_num = 14
doped_base_num = 14

doped_emitter_con = np.power(10, np.linspace(14, 28, doped_emitter_num))
doped_base_con = np.power(10, np.linspace(14, 28, doped_base_num))

X, Y = np.meshgrid(doped_emitter_con, doped_base_con)
eff = pmpp_np / light_source.power_density * 100
eff = ma.masked_where(eff <= 0, eff)
fig2, axes = plt.subplots(2, 2, figsize=(11.25, 8))

cs1 = axes[0, 0].contourf(X, Y, eff, 100, cmap=cm.jet)
axes[0, 0].set_xlabel("Base doping concentration")
axes[0, 0].set_ylabel("Emitter doping concentration")
axes[0, 0].set_yscale("log")
axes[0, 0].set_xscale("log")
axes[0, 0].set_title("Efficiency")

cbar1 = fig2.colorbar(cs1)

cs2 = axes[0, 1].contourf(X, Y, abs(isc_np), 100, cmap=cm.jet)
axes[0, 1].set_xlabel("Base doping concentration")
axes[0, 1].set_ylabel("Emitter doping concentration")
axes[0, 1].set_yscale("log")
axes[0, 1].set_xscale("log")
axes[0, 1].set_title("short circuit current")
cbar2 = fig2.colorbar(cs2)


cs3 = axes[1, 0].contourf(X, Y, abs(voc_np), 100, cmap=cm.jet)
axes[1, 0].set_xlabel("Base doping concentration")
axes[1, 0].set_ylabel("Emitter doping concentration")
axes[1, 0].set_yscale("log")
axes[1, 0].set_xscale("log")
axes[1, 0].set_title("Open circuit voltage")
cbar3 = fig2.colorbar(cs3)

cs4 = axes[1, 1].contourf(X, Y, abs(FF_np) * 100, 100, cmap=cm.jet)
axes[1, 1].set_xlabel("Base doping concentration")
axes[1, 1].set_ylabel("Emitter doping concentration")
axes[1, 1].set_yscale("log")
axes[1, 1].set_xscale("log")
axes[1, 1].set_title("fill factor")
cbar4 = fig2.colorbar(cs4)
plt.show()
