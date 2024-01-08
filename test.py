import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from solcore.light_source import LightSource
from solcore import material
from solcore.constants import electron_mass
from solcore.quantum_mechanics import kp_bands
from solcore import siUnits as unit
from solcore.structure import Structure
from solcore.solar_cell import SolarCell, default_GaAs
from numpy import ma
from matplotlib import cm, ticker

#definition parameter
E_eV = 1.6*10**(-19)
T = 298
wl = np.linspace(350, 1050, 301) * 1e-9
light_source = LightSource(
    source_type="standard",
    version="AM1.5g",
    x=wl,
    output_units="photon_flux_per_m",
    concentration=1,
)
def nkplot(Material,range_min_nm,range_max_nm,plot=False,):
    import matplotlib.pyplot as plt
    from solcore import siUnits as unit
    n = [Material.n(unit(i,'nm')) for i in range(range_min_nm, range_max_nm)]
    k = [Material.k(unit(i,'nm')) for i in range(range_min_nm, range_max_nm)]
    if plot is True:
        plt.plot(n, label=F'n_{Material}', color='red')
        plt.plot(k, label=F'k_{Material}', color='blue')
        plt.legend()
        plt.show()
    return n,k,[i for i in range(range_min_nm, range_max_nm)]


with open('data.npy', 'rb') as fin:
    isc_np = np.load(fin)
    voc_np = np.load(fin)
    FF_np = np.load(fin)
    pmpp_np = np.load(fin)

doped_emitter_num = 14
doped_base_num = 14

doped_emitter_con = np.power(10, np.linspace(14, 28, doped_emitter_num))
doped_base_con = np.power(10, np.linspace(14, 28, doped_base_num))
vint = np.linspace(-2, 1, 300)
V = np.linspace(-2, 0, 150)
X, Y = np.meshgrid(doped_emitter_con, doped_base_con)
eff = pmpp_np / light_source.power_density * 100
eff = ma.masked_where(eff <= 0, eff)

def output(plot, condition_max,**kwargs):
    if plot:
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
        print(doped_emitter_con[x])
        print(doped_base_con[y])
        # try:
        #     print(x2,y2)
        #     print(doped_emitter_con[x2])
        #     print(doped_base_con[y2])
        # except:
        #     pass
        # s = [[str(e) for e in row] for row in eff_list]
        # lens = [max(map(len, col)) for col in zip(*s)]
        # fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
        # table = [fmt.format(*row) for row in s]
        # print('\n'.join(table))
output(True, True)
# print(pmpp_np)
# print(np.nanmax(pmpp_np))
# arr = [22,34,88,99]
# arr1 = [35,43,55,78]
# arr2 = np.nanmax(arr)
# print(arr2)
# GaAs = material('GaAs')()
# nkplot(GaAs,100,1800,plot=True)#add lebel
