import numpy as np
import matplotlib.pyplot as plt

from solcore.light_source import LightSource

# task 1.3 Photon energy in maximum of blackbody radiation at Ts
# Calculate the photon energy at which the black body radiation at 5800 K has a maximum
def task_1_3(bb):
    # print(max(bb.spectrum()[1]))
    # print(bb.spectrum()[1][5])
    count = 0
    for i in bb.spectrum()[0]:
        if bb.spectrum()[1][count] >= max(bb.spectrum()[1]) :
            print(i)
        count += 1
# task 1.4 wavelength of light in the maximum of blackbody radiation at Ts
# Calculate the photon energy at which the black body radiation at 5800 K has a maximum
def task_1_4(bb):
    def nm2eV(wavelenght):
        return 1240/wavelenght
    count = 0
    for i in bb.spectrum()[0]:
        if bb.spectrum()[1][count] >= max(bb.spectrum()[1]):
            print(nm2eV(i))
        count += 1

wl = np.linspace(300, 3000, 200)
bb = LightSource(source_type='black body', x=wl, T=5800, entendue='Sun')

plt.figure(1)
plt.plot(*bb.spectrum(), label='Black body')

plt.xlim(300, 3000)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Power density (Wm$^{-2}$nm$^{-1}$)')
plt.tight_layout()
plt.legend()

# print(best)
task_1_3(bb)
task_1_4(bb)
# plt.show()

