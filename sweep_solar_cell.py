from solcore import material
from solcore.structure import Layer, Junction, TunnelJunction
from solcore.solar_cell import SolarCell
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource
import solcore.poisson_drift_diffusion as PDD
import numpy as np
import matplotlib.pyplot as plt
import time
def GaAs(i):
    # start_time = time.time()
    vint = np.linspace(-2, 1, 300)
    V = np.linspace(-2, 0, 150)
    GaAs_junction = Junction(
        [
            # Layer(width=10e-9, material=window_bottom, role="Window"),  #In width Nd
            Layer(width=i, material=n_GaAs, role="Emitter"),       #Nd width
            Layer(width=2000e-9, material=p_GaAs, role="Base"),         #Na width
            # Layer(width=200e-9, material=bsf_bottom, role="BSF"),       #Na In
        ],
        # sn=1e6,
        # sp=1e6,
        T=T,
        kind="PDD")
    my_solar_cell = SolarCell([GaAs_junction], T=T, substrate=n_GaAs, )
    # print(my_solar_cell)
    solar_cell_solver(my_solar_cell, "iv",
                      user_options={
                          "light_source": light_source,
                          "wavelength": wl,
                          "optics_method": 'TMM',
                          "light_iv": True,
                          "mpp": True,
                          "voltages": V,
                          "internal_voltages": vint, }, )
    end_time = time.time()
    # elapsed_time = end_time - start_time
    # print('-------------------------------------------------------------------------------')
    # print('-------------------------------------------------------------------------------')
    # print('-------------------------------------------------------------------------------')
    #
    # print("Elapsed time: ", elapsed_time)
    # print('-------------------------------------------------------------------------------')
    # print('-------------------------------------------------------------------------------')
    # print('-------------------------------------------------------------------------------')

    return my_solar_cell

T = 298
wl = np.linspace(350, 1050, 301) * 1e-9
window_bottom = material("GaInP")(T=T, Nd=5e24, In=0.49)
n_GaAs = material("GaAs")(T=T, Nd=1e24)
p_GaAs = material("GaAs")(T=T, Na=8e22)
bsf_bottom = material("GaInP")(T=T, Na=5e24, In=0.49)

light_source = LightSource(source_type="standard",
                           version="AM1.5g",x=wl,
                           output_units="photon_flux_per_m",
                           concentration=1)

num_con = 20
con = np.power(10, np.linspace(-9, -4, num_con))
# print(con)
vint = np.linspace(-2, 1, 300)
V = np.linspace(-2, 0, 150)
# for i in con:
#     print(f'{i:e}')

allI = []
isc = []
voc = []
FF = []
pmpp = []
fig3, axIV = plt.subplots(1, 1, figsize=(6, 4))
cm = plt.get_cmap('gist_rainbow')

# my_solar_cell = GaAs(150e-9)
# isc.append(my_solar_cell.iv["Isc"])
# voc.append(my_solar_cell.iv["Voc"])
# FF.append(my_solar_cell.iv["FF"])
# pmpp.append(my_solar_cell.iv["Pmpp"])
# allI.append(my_solar_cell.iv["IV"][1])

def seq():
    my_solar_cell = GaAs(1e-9)
    isc.append(my_solar_cell.iv["Isc"])
    voc.append(my_solar_cell.iv["Voc"])
    FF.append(my_solar_cell.iv["FF"])
    pmpp.append(my_solar_cell.iv["Pmpp"])
    allI.append(my_solar_cell.iv["IV"][1])
    axIV.plot(-V, my_solar_cell.iv["IV"][1] / -1e2, label=f'emitter width{1e-9}')
    my_solar_cell = GaAs(10e-9)
    isc.append(my_solar_cell.iv["Isc"])
    voc.append(my_solar_cell.iv["Voc"])
    FF.append(my_solar_cell.iv["FF"])
    pmpp.append(my_solar_cell.iv["Pmpp"])
    allI.append(my_solar_cell.iv["IV"][1])
    axIV.plot(-V, my_solar_cell.iv["IV"][1] / -1e2, label=f'emitter width{10e-9}')
    my_solar_cell = GaAs(100e-9)
    isc.append(my_solar_cell.iv["Isc"])
    voc.append(my_solar_cell.iv["Voc"])
    FF.append(my_solar_cell.iv["FF"])
    pmpp.append(my_solar_cell.iv["Pmpp"])
    allI.append(my_solar_cell.iv["IV"][1])
    axIV.plot(-V, my_solar_cell.iv["IV"][1] / -1e2, label=f'emitter width{100e-9}')
    my_solar_cell = GaAs(1000e-9)
    isc.append(my_solar_cell.iv["Isc"])
    voc.append(my_solar_cell.iv["Voc"])
    FF.append(my_solar_cell.iv["FF"])
    pmpp.append(my_solar_cell.iv["Pmpp"])
    allI.append(my_solar_cell.iv["IV"][1])
    axIV.plot(-V, my_solar_cell.iv["IV"][1] / -1e2, label=f'emitter width{1000e-9}')
    axIV.legend(loc="lower left", frameon=False)
    axIV.set_ylim(10, 0)
    axIV.set_xlim(0, 1.5)
    axIV.set_xlabel("Voltage (V)")
    axIV.set_ylabel("amp/cm^2")
    plt.show()
def loop():
    axIV.set_prop_cycle(color=[cm(1. * i / num_con) for i in range(num_con)])
    start_time = time.time()
    for i in con:
        my_solar_cell = GaAs(i)
        isc.append(my_solar_cell.iv["Isc"])
        voc.append(my_solar_cell.iv["Voc"])
        FF.append(my_solar_cell.iv["FF"])
        pmpp.append(my_solar_cell.iv["Pmpp"])
        allI.append(my_solar_cell.iv["IV"][1])
        # And now, everything is plotting...
        axIV.plot(-V, my_solar_cell.iv["IV"][1] / -1e2, label=f'emitter width{i:e}')
    axIV.legend(loc="lower right", frameon=False)
    axIV.set_ylim(3, 0)
    axIV.set_xlim(0, 1.5)
    axIV.set_xlabel("Voltage (V)")
    axIV.set_ylabel("amp/cm^2")
    plt.tight_layout()
    end_time = time.time()
    elapsed_time = end_time - start_time
    print('-------------------------------------------------------------------------------')
    print('-------------------------------------------------------------------------------')
    print('-------------------------------------------------------------------------------')

    print("Elapsed time: ", elapsed_time)
    print('-------------------------------------------------------------------------------')
    print('-------------------------------------------------------------------------------')
    print('-------------------------------------------------------------------------------')
    plt.show()


# seq()
loop()