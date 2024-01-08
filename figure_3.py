from solcore import material
from solcore.structure import SolcoreMaterialToStr
from solcore.quantum_mechanics import kp_bands
from solcore.material_system.material_system import nkdb_load_n
import matplotlib.pyplot as plt

GaAs = material('GaAs')()
GaAs_sopra = material("GaAs", sopra=True)()
Ge   = material('Ge')()
def nkplot(Material,range_min_nm,range_max_nm,plot=True,linestyle='solid'):
    import matplotlib.pyplot as plt
    from solcore import siUnits as unit
    n = [Material.n(unit(i,'nm')) for i in range(range_min_nm, range_max_nm)]
    k = [Material.k(unit(i,'nm')) for i in range(range_min_nm, range_max_nm)]
    if plot is True:
        plt.plot(n, label=F'n_{Material}', color='red',linestyle=linestyle)
        plt.plot(k, label=F'k_{Material}', color='blue',linestyle=linestyle)
        plt.legend()
        # plt.show()
    return n, k, [i for i in range(range_min_nm, range_max_nm)]
# def absorb(mat):
#     import
nkplot(GaAs,300,1800)
nkplot(Ge,300,1800,linestyle='dashed')
plt.legend(loc='upper right')
print(SolcoreMaterialToStr(GaAs))
print(SolcoreMaterialToStr(Ge))
plt.tight_layout()
plt.show()