import solcore
from solcore import get_parameter
from solcore import si, material
import solcore.quantum_mechanics as QM

from solcore.structure import Layer, Structure, Junction

# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    T = 150
    n_AlGaAs = material('AlGaAs')(T=T,Al=0.3,Nd=1e25)
    p_AlGaAs = material('AlGaAs')(T=T,Al=0.3,Na=1e25)
    i_GaAs = material('GaAs')(T=T)
    struc = Structure([
        Layer(width=1e-6, material=n_AlGaAs),
        Layer(width=0.1e-6, material=i_GaAs),
        Layer(width=1e-6, material=p_AlGaAs)
    ], substrate = i_GaAs, T=T)
    SR, band= QM.schrodinger(struc, quasiconfined=0.05, graphtype='potentialsLDOS', num_eigenvalues=200,show=True,periodic=False)


# See PyCharm help at https://www.jetbrains.com/help/pycharm/
