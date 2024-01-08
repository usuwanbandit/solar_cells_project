from solcore import si, material
from solcore.structure import Layer, Structure, Junction
from solcore.solar_cell import SolarCell
import solcore.quantum_mechanics as QM
import solcore.poisson_drift_diffusion as PDD
import numpy as np
from solcore.light_source import LightSource
from solcore.solar_cell_solver import solar_cell_solver
import matplotlib.pyplot as plt
import os
import shutil
def create_folder(create_folder):
    import os
    current_path = os.getcwd()
    # print(current_path)
    current_path = os.path.join(current_path, create_folder)
    print(current_path)
    if not os.path.exists(create_folder):
        os.makedirs(create_folder)
        print('create folder success')


def save_tuple2text(text, name, *data):  # each data ==> tuple(name, detail)
    text += name + '\nCustum text version' + '\n'
    text += "===============================START======================================= \n"
    for i in data:
        text += str(i) +'\n'
    text += "================================END======================================= \n"
    return text  # note data topic and keyword


def save_full_solar_cells(text, solar_cells, name_solar):
    text += name_solar + 'Full text version' + '\n'
    text += "===============================START======================================= \n"
    for i in solar_cells:
        text += str(i) + '\n'
    text += "================================END======================================= \n"
    return text


def save_file_direction(save_folder, text, name_text):  # find from current file
    import os
    current_path = os.getcwd()
    current_path = os.path.join(current_path, save_folder)
    if not os.path.exists(current_path):
        os.makedirs(current_path)
        print(f'create {current_path} folder ')
    complete_Name = os.path.join(current_path, name_text + ".txt")
    file1 = open(complete_Name, "w")
    toFile = text
    file1.write(toFile)
    file1.close()
    print('save success')


def showSLQW():
    T = 293
    n_GaAs = material("GaAs")(T=T, Nd=1e24)
    p_GaAs = material("GaAs")(T=T, Na=8e22)
    barrier = material('AlGaAs')(T=T, Al=0.3, Nd=1e18)
    QWmat1 = material("InSb")(T=T, strained=True)
    QWmat2 = material("GaSb")(T=T, strained=True)
    i_GaAs = material("GaAs")(T=T)
    struc = Structure([
                            Layer(width=100e-9, material=barrier, role='barrier'),
                            Layer(width=83.84e-9, material=i_GaAs, role="interlayer"),
                            Layer(width=16.16e-9, material=QWmat1, role="well"),
                            Layer(width=96.46e-9, material=i_GaAs, role="interlayer"),
                            Layer(width=2.54e-9, material=QWmat2, role="well"),
                            Layer(width=50e-9, material=n_GaAs, role="interlayer"),
                            Layer(width=100e-9, material=barrier, role="barrier"),
                        ],substrate=p_GaAs, T=T)
    # print(struc)
    SR, band= QM.schrodinger(struc, quasiconfined=0, graphtype='potentialsLDOS', num_eigenvalues=200,show=True,periodic=False)
    # print(output_2)
    return SR, band
# text = ''
# solar_text = save_full_solar_cells('', my_solar_cell, 'test_solar')
# save_file_direction('test_folder', solar_text, 'test1')
wl = np.linspace(350, 1200, 401) * 1e-9
light_source = LightSource(
    source_type="standard",
    version="AM1.5g",
    x=wl,
    output_units="photon_flux_per_m",
    concentration=1,
)

def showeiei(barrier, interlayer, dot):
    T = 293
    n_GaAs = material("GaAs")(T=T, Nd=1e24)
    p_GaAs = material("GaAs")(T=T, Na=8e22)
    #dot type 1
    # QWmat = material("InGaAs")(T=T, In=0.2, strained=True)
    # Bmat = material("GaAsP")(T=T, P=0.1, strained=False)
    # i_GaAs = material("GaAs")(T=T)
    #dottype2
    Bmat = material('AlGaAs')(T=T, Al=0.3, Nd=1e18)
    QWmat = material("InSb")(T=T, strained=True)
    QWmat2 = material("GaSb")(T=T, strained=True)
    i_GaAs = material("GaAs")(T=T)
    # i_GaAs_dope = material("GaAs")(T=T, Na= 1e14)

    struc = Structure(
            # [
                # Layer(width=50e-9, material=n_GaAs, role="Emitter")]+
                      # [Layer(width=barrier, material=Bmat, role="barrier")] +
                      #dot type 1---------------------------------------
                        # [Layer(width=barrier, material=Bmat, role="barrier"),
                        #       Layer(width=interlayer, material=i_GaAs, role="interlayer"),
                        #       Layer(width=dot, material=QWmat, role="well"),
                        #       Layer(width=interlayer, material=i_GaAs, role="interlayer"),
                        #       Layer(width=barrier, material=Bmat, role="barrier")]
                            #   +
                        #dot type 2
                         [Layer(width=barrier, material=Bmat, role="barrier"),
                              Layer(width=interlayer, material=i_GaAs, role="interlayer"),
                              Layer(width=dot, material=QWmat, role="well"),
                              Layer(width=interlayer, material=i_GaAs, role="interlayer"),
                              Layer(width=barrier, material=Bmat, role="barrier")]
                      # +10 * [Layer(width=10e-9, material=Bmat, role="barrier"),
                             # Layer(width=2e-9, material=i_GaAs, role="well"),
                             # Layer(width=7e-9, material=QWmat, role="well"),
                             # Layer(width=2e-9, material=i_GaAs, role="well"),
                             # Layer(width=10e-9, material=Bmat, role="barrier"),] +
                    #   [Layer(width=barrier, material=Bmat, role="barrier")] +
                    #   [Layer(width=50e-9, material=p_GaAs, role="Base")
                    #    ]
                      ,substrate=i_GaAs)
    # print(struc)
    SR, band = QM.schrodinger(struc, quasiconfined=0, graphtype='potentialsLDOS', num_eigenvalues=200, show=True)
    # print(output_2)
    return SR, band
def plotice(out, band):
    plt.contourf(out['LDOS']['x'] * 1e9, out['LDOS']['Ee'] / 1e-19, out['LDOS']['LDOSe'], 100, cmap='gnuplot2_r', vmin=0, vmax=max(out['LDOS']['LDOSe'].flatten()) * 1.2)
    plt.plot(band['x']*1e9, band['Ve']/1e-19, 'k', label="Ve")
    plt.plot(band['x']*1e9, band['Vlh']/1e-19, 'k--',  label="Vlh")
    plt.plot(band['x']*1e9, band['Vhh']/1e-19, 'k',  label="Vhh")

    plt.contourf(out['LDOS']['x'] * 1e9, out['LDOS']['Eh'] / 1e-19, out['LDOS']['LDOSh'], 100, cmap='gnuplot2_r', vmin=0, vmax=max(out['LDOS']['LDOSh'].flatten()) * 1.2)
    plt.ylabel('Energy (eV)')
    plt.xlabel('Possition (nm)')
    plt.legend()
    plt.show()

def solar_cells(barrier, interlayer, dot, update, version,active_dot=True):

    T = 298
    wl = np.linspace(350, 1050, 301) * 1e-9

    # First, we create the materials of the QW
    # สร้างวัสดุทำรับQW ในsolar cell
    QWmat = material("InGaAs")(T=T, In=0.2, strained=True)
    Bmat = material("GaAsP")(T=T, P=0.1, strained=False)
    i_GaAs = material("GaAs")(T=T)
    i_GaAs_dope = material("GaAs")(T=T, Na=1e14)

    n_GaAs = material("GaAs")(T=T, Nd=1e24)
    p_GaAs = material("GaAs")(T=T, Na=8e22)
    if active_dot:
        QW = PDD.QWunit([
            Layer(width=barrier, material=Bmat, role="barrier"),
            Layer(width=interlayer, material=i_GaAs, role="interlayer"),
            Layer(width=dot, material=QWmat, role="well"),
            Layer(width=interlayer, material=i_GaAs, role="interlayer"),
            Layer(width=barrier, material=Bmat, role="barrier")
        ], T=T, repeat=30, substrate=i_GaAs, )
        QW_list = QW.GetEffectiveQW(wavelengths=wl, mode="kp6x6")
        GaAs_junction = Junction(
            [Layer(width=150e-9, material=n_GaAs, role="Emitter"),]
             # Layer(width=barrier, material=Bmat, role="barrier")]
            + QW_list
            # + [Layer(width=barrier, material=Bmat, role="barrier"),
            +[Layer(width=2000e-9, material=p_GaAs, role="Base"),
            ],
            sn=1e6,
            sp=1e6,T=T,kind="PDD",)
    else:
        GaAs_junction = Junction(
            [Layer(width=150e-9, material=n_GaAs, role="Emitter"), ]
            # Layer(width=barrier, material=Bmat, role="barrier")]
            # + QW_list
                # + [Layer(width=barrier, material=Bmat, role="barrier"),
            +[Layer(width=2000e-9, material=p_GaAs, role="Base"),
            ],
            sn=1e6,
            sp=1e6, T=T, kind="PDD", )
    MgF2 = material("MgF2")()
    ZnS = material("ZnScub")()
    my_solar_cell = SolarCell([
        # [Layer(width=110e-9, material=MgF2, role="ARC1"),
        # Layer(width=60e-9, material=ZnS, role="ARC2")
                            GaAs_junction,]
                              ,T=T,substrate=p_GaAs,)
    solar_cell_solver(my_solar_cell, "qe",
                      user_options={"light_source": light_source,
                                    "wavelength": wl,
                                    "optics_method": "TMM",}, )

    num_con = 5  # จำนวนในการสร้างแสง
    con = np.logspace(0, 3, num_con)
    vint = np.linspace(-3.5, 4, 600)
    V = np.linspace(-3.5, 0, 300)
    allI = []; isc = []; voc = []; FF = []; pmpp = []
    fig3, axIV = plt.subplots(1, 1, figsize=(6, 4))
    for i,c in enumerate(con):  # ทำการยิงแสงทั้งหมด20ครั้งตามcon
        light_source.concentration = c
        solar_cell_solver(my_solar_cell,"iv"
                          ,user_options={"light_source": light_source,
                                        "wavelength": wl,
                                         "optics_method": None,
                                         "light_iv": True,
                                         "mpp": True,
                                         "voltages": V,
                                         "internal_voltages": vint,
                            },)
        isc.append(my_solar_cell.iv["Isc"])
        voc.append(my_solar_cell.iv["Voc"])
        FF.append(my_solar_cell.iv["FF"])
        pmpp.append(my_solar_cell.iv["Pmpp"])
        allI.append(my_solar_cell.iv["IV"][1])
        axIV.plot(-V, my_solar_cell.iv["IV"][1] / isc[-1], label=int(c))  # เป็นการวัดลักษณะIVเมิอแสงเปลี่ยนไป
        if active_dot:
            print('======================================================================')
            print('======================================================================')
            print(f'{(i+num_con*update[0])/(num_con*update[1])*100}%')
            print('======================================================================')
            print('======================================================================')
        else:
            print('======================================================================')
            print('======================================================================')

            print(f'{i}/{num_con} test_case ')
            print('======================================================================')
            print('======================================================================')

    with open('data.npy', 'wb') as fout:
        np.save(fout, np.array(isc))
        np.save(fout, np.array(voc))
        np.save(fout, np.array(FF))
        np.save(fout, np.array(pmpp))
        np.save(fout, con)
    axIV.legend(loc="lower left", frameon=False)
    axIV.set_ylim(0, 1.1)
    axIV.set_xlim(0, 1.5)
    axIV.set_xlabel("Voltage (V)")
    axIV.set_ylabel("Normalized current (-)")
    plt.tight_layout()

    fig2, axes = plt.subplots(2, 2, figsize=(11.25, 8))

    axes[0, 0].semilogx(con, np.array(pmpp) / con / 10,"r-o")  # ทำไมต้องหาร10ด้วยไม่เข้าใจเลย(ไม่ใช้หารด้วย 0.97หรอเพราะว่าAM1.5g = 970w/m^2)แล้วทำไมไม่x100
    axes[0, 0].set_xlabel("Concentration (suns)")
    axes[0, 0].set_ylabel("Efficiency (%)")

    axes[0, 1].loglog(con, abs(np.array(isc)), "b-o")
    axes[0, 1].set_xlabel("Concentration (suns)")
    axes[0, 1].set_ylabel("I$_{SC}$ (Am$^{-2}$)")

    axes[1, 0].semilogx(con, abs(np.array(voc)), "g-o")
    axes[1, 0].set_xlabel("Concentration (suns)")
    axes[1, 0].set_ylabel("V$_{OC}$ (V)")

    axes[1, 1].semilogx(con, abs(np.array(FF)) * 100, "k-o")
    axes[1, 1].set_xlabel("Concentration (suns)")
    axes[1, 1].set_ylabel("Fill Factor (%)")

    plt.tight_layout()
    fig, ax1  = plt.subplots(1, 1, figsize=(6, 4))
    # We can plot the electron and hole densities in equilibrium and at short circuit
    ax1.plot(wl * 1e9, my_solar_cell.absorbed, "k", label="Total Absorbed")
    ax1.legend(loc="upper right", frameon=False)
    ax1.set_xlabel("Wavelength (nm)")
    ax1.set_ylabel("EQE")
    ax1.set_ylim(0, 1.1)
    ax1.set_xlim(350, 1150)
    plt.tight_layout()

    fig3.savefig(f'IV_curve_QW_{dot:.2e}_width_{version}.png', dpi=300)
    fig2.savefig(f'performance_QW_{dot:.2e}_width_{version}.png', dpi=300)
    fig.savefig(f'EQE_QW_{dot:.2e}_width_{version}.png', dpi=300)
    text = ''
    data = save_tuple2text(text,'profile'+str(dot)+'width', my_solar_cell)
    for x in my_solar_cell:
        data += str(x)
    save_file_direction(f'data_of_dot_width{dot:.2e}_{version}', data, f'QWcapping_width{dot:.2e}_{version}')
    def movefile(file,direction):
        save_path = os.path.join(current_path, direction)
        fig1_loc = os.path.join(current_path, file)
        fig1_loc_new = os.path.join(save_path, file)
        shutil.move(fig1_loc, fig1_loc_new)

    current_path = os.getcwd()
    movefile(f'EQE_QW_{dot:.2e}_width_{version}.png', f'data_of_dot_width{dot:.2e}_{version}')
    movefile(f'performance_QW_{dot:.2e}_width_{version}.png', f'data_of_dot_width{dot:.2e}_{version}')
    movefile(f'IV_curve_QW_{dot:.2e}_width_{version}.png', f'data_of_dot_width{dot:.2e}_{version}')

def loadgraph(filename):
    with open(filename, 'rb') as fin:
        isc = np.load(fin)
        voc = np.load(fin)
        FF = np.load(fin)
        pmpp = np.load(fin)
        con = np.load(fin)
    fig2, axes = plt.subplots(2, 2, figsize=(11.25, 8))

    axes[0, 0].semilogx(con, np.array(pmpp) / con / 10,
                        "r-o")  # ทำไมต้องหาร10ด้วยไม่เข้าใจเลย(ไม่ใช้หารด้วย 0.97หรอเพราะว่าAM1.5g = 970w/m^2)แล้วทำไมไม่x100
    axes[0, 0].set_xlabel("Concentration (suns)")
    axes[0, 0].set_ylabel("Efficiency (%)")

    axes[0, 1].loglog(con, abs(np.array(isc)), "b-o")
    axes[0, 1].set_xlabel("Concentration (suns)")
    axes[0, 1].set_ylabel("I$_{SC}$ (Am$^{-2}$)")

    axes[1, 0].semilogx(con, abs(np.array(voc)), "g-o")
    axes[1, 0].set_xlabel("Concentration (suns)")
    axes[1, 0].set_ylabel("V$_{OC}$ (V)")

    axes[1, 1].semilogx(con, abs(np.array(FF)) * 100, "k-o")
    axes[1, 1].set_xlabel("Concentration (suns)")
    axes[1, 1].set_ylabel("Fill Factor (%)")

    plt.tight_layout()
    plt.show()

def eqe(barrier, interlayer, dot):
    T = 298
    wl = np.linspace(350, 1050, 301) * 1e-9

    # First, we create the materials of the QW
    # สร้างวัสดุทำรับQW ในsolar cell
    QWmat = material("InGaAs")(T=T, In=0.2, strained=True)
    Bmat = material("GaAsP")(T=T, P=0.1, strained=True)
    i_GaAs = material("GaAs")(T=T)
    i_GaAs_dope = material("GaAs")(T=T, Na=1e14)

    n_GaAs = material("GaAs")(T=T, Nd=1e24)
    p_GaAs = material("GaAs")(T=T, Na=8e22)
    QW = PDD.QWunit([
        Layer(width=barrier, material=Bmat, role="barrier"),
        Layer(width=interlayer, material=i_GaAs, role="interlayer"),
        Layer(width=dot, material=QWmat, role="well"),
        Layer(width=interlayer, material=i_GaAs, role="interlayer"),
        Layer(width=barrier, material=Bmat, role="barrier")
    ], T=T, repeat=30, substrate=i_GaAs, )
    QW_list = QW.GetEffectiveQW(wavelengths=wl)
    GaAs_junction = Junction(
        [Layer(width=150e-9, material=n_GaAs, role="Emitter"),]
         # Layer(width=barrier, material=Bmat, role="barrier")]
        + QW_list
        # + [Layer(width=barrier, material=Bmat, role="barrier"),
        +[Layer(width=2000e-9, material=p_GaAs, role="Base"),
        ],
        sn=1e6,
        sp=1e6,T=T,kind="PDD",)
    MgF2 = material("MgF2")()
    ZnS = material("ZnScub")()
    my_solar_cell = SolarCell([Layer(width=110e-9, material=MgF2, role="ARC1"),
                               Layer(width=60e-9, material=ZnS, role="ARC2")
                            ,GaAs_junction,]
                              ,T=T,substrate=n_GaAs,)
    solar_cell_solver(my_solar_cell, "qe",
                      user_options={"light_source": light_source,
                                    "wavelength": wl,
                                    "optics_method": "TMM",}, )
    GaAs_junction1 = Junction(
        [Layer(width=150e-9, material=n_GaAs, role="Emitter"),]
         # Layer(width=barrier, material=Bmat, role="barrier")]
        # + QW_list
        # + [Layer(width=barrier, material=Bmat, role="barrier"),
        +[Layer(width=2000e-9, material=p_GaAs, role="Base"),
        ],
        sn=1e6,
        sp=1e6,T=T,kind="PDD",)
    MgF2 = material("MgF2")()
    ZnS = material("ZnScub")()
    my_solar_cell1 = SolarCell([Layer(width=110e-9, material=MgF2, role="ARC1"),
                               Layer(width=60e-9, material=ZnS, role="ARC2")
                            ,GaAs_junction1,]
                              ,T=T,substrate=n_GaAs,)
    solar_cell_solver(my_solar_cell1, "qe",
                      user_options={"light_source": light_source,
                                    "wavelength": wl,
                                    "optics_method": "TMM",}, )
    plt.plot(wl * 1e9, my_solar_cell.absorbed * 100, "k", label="GaAs_with_QD", color='r')
    plt.plot(wl * 1e9, my_solar_cell1.absorbed * 100, "k", label="GaAs", color='b')
    plt.legend(loc="upper right", frameon=False)
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("EQE (%)")
    plt.ylim(0, 110)
    plt.xlim(350, 1150)
    plt.show()

# showeiei2()
# SR, band = showSLQW()
# plotice(SR, band)
# loadgraph('data.npy')
# solar_cells(5e-9, 2e-9,7e-9, (1,1),'powerpoint_show', active_dot=True)
#--------------------------------------------
showSLQW()
#--------------------------------------------

#---------------------------------
# showeiei(100e-9,100e-9,16.16e-9)
# showeiei(100e-9,1000e-9,60e-9)
#--------------------------------------------

# con_num = 8
# con = np.linspace(3e-9, 10e-9,con_num)
# eqe(4e-9, 3e-9, 7e-9)
# for j,i in enumerate(con):
#     # print()
#     solar_cells(5e-9, 2e-9,i, (j,con_num),'capping5nm_ver0.2', active_dot=True)
#     # showeiei(2e-9,2e-9,i)
#     print('======================================================================')
#     print(f'{j+1}/{con_num}', 'save susses')
#     print('======================================================================')
# solar_cells()
