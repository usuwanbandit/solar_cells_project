import solcore
from solcore import get_parameter
from solcore import si, material
import solcore.quantum_mechanics as QM
import numpy as np
from solcore.structure import Layer, Structure, Junction

# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.



# Press the green button in the gutter to run the script.
# if __name__ == '__main__':
#     T = 150
#     n_AlGaAs = material('AlGaAs')(T=T,Al=0.3,Nd=1e25)
#     p_AlGaAs = material('AlGaAs')(T=T,Al=0.3,Na=1e25)
#     i_GaAs = material('GaAs')(T=T)
#     struc = Structure([
#         Layer(width=1e-6, material=n_AlGaAs),
#         Layer(width=0.1e-6, material=i_GaAs),
#         Layer(width=1e-6, material=p_AlGaAs)
#     ], substrate = i_GaAs, T = T)
#     SR, band= QM.schrodinger(struc, quasiconfined=0.05, graphtype='potentialsLDOS', num_eigenvalues=200,show=True,periodic=False)
T=150
def linecomp(start,stop,num,width):
    list_run = np.linspace(start,stop,num)/100
    list_layer = [Layer(width=width/num, material=material('AlGaAs')(T=T, Al=i)) for i in list_run]
    print(list_layer)
    return list_layer
def laser():
    T = 150
    n_AlGaAs = material('AlGaAs')(T=T,Al=0.45,Nd=1e25)
    p_AlGaAs = material('AlGaAs')(T=T,Al=0.45,Na=1e25)
    n_laAlGaAs = linecomp(45,20,1000,120e-9)
    p_laAlGaAs = linecomp(20,45,1000,120e-9)
    i_GaAs = material('GaAs')(T=T)
    struc = Structure([
        Layer(width=1e-6, material=n_AlGaAs)]
        +n_laAlGaAs
        +[Layer(width=10e-9, material=i_GaAs)]
        +p_laAlGaAs
        +[Layer(width=0.1e-6, material=p_AlGaAs),
        Layer(width=1e-6, material=p_AlGaAs)
    ], substrate = i_GaAs, T = T)
    SR, band= QM.schrodinger(struc, quasiconfined=0.05, graphtype='potentialsLDOS', num_eigenvalues=200,show=True,periodic=False)
#
laser()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
