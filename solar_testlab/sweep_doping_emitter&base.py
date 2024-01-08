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
light_source = LightSource(
    source_type="standard",
    version="AM1.5g",
    x=wl,
    output_units="photon_flux_per_m",
    concentration=1,
)
def solver_GaAs_for_sweep(Nd, width):
    n_GaAs = material("GaAs")(T=T, Nd=Nd)
    p_GaAs = material("GaAs")(T=T, Na=8e22)
    GaAs_junction = Junction(
    [
        # Layer(width=10e-9, material=window_bottom, role="Window"),  #In width Nd
        Layer(width=width, material=n_GaAs, role="Emitter"),       #Nd width
        Layer(width=200e-9, material=p_GaAs, role="Base"),         #Na width
        # Layer(width=200e-9, material=bsf_bottom, role="BSF"),       #Na In
    ],
    T=T,
    kind="PDD")
    my_solar_cell = SolarCell([
        GaAs_junction,
        Layer(width=20000e-9, material=p_GaAs, role="Base"),  # Na width
    ],T=T,substrate=p_GaAs,)
    solar_cell_solver(my_solar_cell,"iv",
               user_options={
                   "light_source": light_source,
                   "wavelength": wl,
                   "optics_method": 'TMM',
                   "light_iv": True,
                   "mpp": True,
                   "voltages": V,
                   "internal_voltages": vint,
    },)
    return my_solar_cell

vint = np.linspace(-2, 1, 300)
V = np.linspace(-2, 0, 150)
def sweep_width():
    width_con = 14
    # doped_con = 19

    con_width = np.power(10, np.linspace(-9, -4, width_con))
    # con_doped = np.power(10, np.linspace(14, 28, doped_con))
    vint = np.linspace(-2, 1, 300)
    V = np.linspace(-2, 0, 150)

    isc_np = np.zeros(width_con)
    voc_np = np.zeros(width_con)
    FF_np = np.zeros(width_con)
    pmpp_np = np.zeros(width_con)
    # allI_np = np.zero((width_con, width_con))
    index = 0
    N = width_con
    for i, width in enumerate(con_width):
        my_solar_cell = solver_GaAs_for_sweep(1e24, width)
        isc_np[i] = my_solar_cell.iv["Isc"]
        voc_np[i] = my_solar_cell.iv["Voc"]
        FF_np[i] = my_solar_cell.iv["FF"]
        pmpp_np[i] = my_solar_cell.iv["Pmpp"]
        # allI_np[i, j] =my_solar_cell.iv["IV"][1]
        # emitter_width(width)
        index += 1
        print('==============================\n')
        print(int(index / N * 100), "%\n")
        print('==============================\n')
    with open('../data.npy', 'wb') as fout:
        np.save(fout, isc_np)
        np.save(fout, voc_np)
        np.save(fout, FF_np)
        np.save(fout, pmpp_np)
        np.save(fout, con_width)
    # return
def output(plot, condition_max,**kwargs):
    with open('../data.npy', 'rb') as fin:
        isc_np = np.load(fin)
        voc_np = np.load(fin)
        FF_np = np.load(fin)
        pmpp_np = np.load(fin)
        con_width = np.load(fin)
    eff = pmpp_np/light_source.concentration *100
    if plot:
        fig2, axes = plt.subplots(2, 2, figsize=(11.25, 8))
        axes[0, 0].plot(con_width, eff)
        axes[0, 0].set_xlabel("emitter width")
        axes[0, 0].set_ylabel("eff")
        # axes[0, 0].set_yscale("log")
        axes[0, 0].set_xscale("log")
        axes[0, 0].set_title("Efficiency")

        # cbar1 = fig2.colorbar(cs1)

        cs2 = axes[0, 1].plot(con_width, abs(isc_np),)
        axes[0, 1].set_xlabel("emitter width")
        axes[0, 1].set_ylabel("I$_{SC}$ (Am$^{-2}$)")
        # axes[0, 1].set_yscale("log")
        axes[0, 1].set_xscale("log")
        axes[0, 1].set_title("short circuit current")
        # cbar2 = fig2.colorbar(cs2)


        cs3 = axes[1, 0].plot(con_width, abs(voc_np),)
        axes[1, 0].set_xlabel("emitter width")
        axes[1, 0].set_ylabel("V$_{OC}$ (V)")
        # axes[1, 0].set_yscale("log")
        axes[1, 0].set_xscale("log")
        axes[1, 0].set_title("Open circuit voltage")
        # cbar3 = fig2.colorbar(cs3)

        cs4 = axes[1, 1].plot(con_width, abs(FF_np) * 100,)
        axes[1, 1].set_xlabel("emitter width")
        axes[1, 1].set_ylabel("Fill Factor (%)")
        # axes[1, 1].set_yscale("log")
        axes[1, 1].set_xscale("log")
        axes[1, 1].set_title("fill factor")
        # cbar4 = fig2.colorbar(cs4)

        plt.show()
    if condition_max:
        eff_list = np.ndarray.tolist(eff)
        max_eff = 0
        for i, pmpp_y in enumerate(eff_list):
            for j, pmpp_x in enumerate(pmpp_y):

                if max_eff < pmpp_x:
                    print(max_eff)
                    max_eff = pmpp_x
                    x = j
                    y = i
                elif max_eff == pmpp_x:
                    x2 = j
                    y2 = i
        print(max_eff, 'maximum efficiency')
        # print(doped_emitter_con[x])
        # print(doped_base_con[y])
# sweep_width()
# output(True,False)

my_solar_cell = solver_GaAs_for_sweep(8e23, 15e-9)
# print(my_solar_cell)
isc_np = my_solar_cell.iv["Isc"]
voc_np = my_solar_cell.iv["Voc"]
FF_np = my_solar_cell.iv["FF"]
pmpp_np = my_solar_cell.iv["Pmpp"]
print('isc', isc_np)
print('voc', voc_np)
print('FF', FF_np)
print('pmpp', pmpp_np)
# best_solar = solver_GaAs_for_sweep(4.923882631706752e23,4.923882631706752e23)
# print(best_solar.__dict__)
# print(best_solar)
