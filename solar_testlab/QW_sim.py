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
import pickle
def create_folder(create_folder):
    import os
    current_path = os.getcwd()
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

bulk = material("GaAs")(T=293, strained=False)
barrier = material("GaAsP")(T=293, P=0.1, strained=True)

top_layer = Layer(width=si("30nm"), material=barrier)
inter = Layer(width=si("3nm"), material=bulk)
barrier_layer = Layer(width=si("5nm"), material=barrier)
bottom_layer = top_layer
wl = np.linspace(350, 1050, 301) * 1e-9


InSb = material('InSb')(T=293, strained=False)
GaSb = material('GaSb')(T=293, strained=False)
buffer = material('AlGaAs')(T=293, Al=0.3)
capping = material('GaAs')(T=293)
n_GaAs = material('GaAs')(T=293, Nd=1e15)
p_GaAs = material('GaAs')(T=293, Na=1e15)
QW = material("InGaAs")(T=293, In=0.15, strained=True)

n_layer = Layer(width=250e-9, material=n_GaAs)
QW1 = Layer(width=5e-9, material=InSb)
Capping_layer = Layer(width=80e-9, material=capping)
QW2 = Layer(width=5e-9, material=GaSb)
p_layer = Layer(width=400e-9, material=p_GaAs)
test_ver1 = Structure([n_layer, QW1, Capping_layer, QW2, p_layer], substrate=bulk)
test_ver2 = Structure([n_layer, QW1, Capping_layer, p_layer], substrate=bulk)
def old_met():
    # outputtest1 = QM.schrodinger(test_ver1, quasiconfined=0, graphtype='potentials', num_eigenvalues=20,show=True)
    # outputtest2 = QM.schrodinger(test_ver2, quasiconfined=0, graphtype='potentialsLDOS', num_eigenvalues=200,show=True)

    # And the layer
    well_layer = Layer(width=si("7.2nm"), material=QW)

    # The following lines create the QW structure, with different number of QWs and interlayers. Indicating the substrate
    # material with the keyword "substrate" is essential in order to calculate correctly the strain.

    # A single QW with interlayers
    test_structure_1 = Structure([top_layer, inter, well_layer, inter, bottom_layer], substrate=bulk)
    # output_1 = QM.schrodinger(test_structure_1, quasiconfined=0, graphtype='potentials', num_eigenvalues=20, show=True)
    # output_1 = QM.schrodinger(test_structure_1, quasiconfined=0, graphtype='potentialsLDOS', num_eigenvalues=20, show=True)
    # 10 QWs without interlayers
    test_structure_2 = Structure([top_layer, barrier_layer] + 10 * [well_layer, barrier_layer] + [bottom_layer], substrate=bulk)

    output_2 = QM.schrodinger(test_structure_2, quasiconfined=0.05, graphtype='potentialsLDOS', num_eigenvalues=200,show=True)
def showschrodinger():
    T = 293
    n_GaAs = material("GaAs")(T=T, Nd=1e24)
    p_GaAs = material("GaAs")(T=T, Na=8e22)
    QWmat = material("InGaAs")(T=T, In=0.2, strained=True)
    Bmat = material("GaAsP")(T=T, P=0.1, strained=True)
    i_GaAs = material("GaAs")(T=T)
    i_GaAs_dope = material("GaAs")(T=T, Na= 1e14)

    struc = Structure([Layer(width=50e-9, material=n_GaAs, role="Emitter")]+
                      [Layer(width=2e-9, material=Bmat, role="barrier")] +
                      + 10 * [Layer(width=2e-9, material=Bmat, role="barrier"),
                              Layer(width=1e-9, material=i_GaAs_dope, role="well"),
                              Layer(width=7e-9, material=QWmat, role="well"),
                              Layer(width=1e-9, material=i_GaAs_dope, role="well"),
                              Layer(width=2e-9, material=Bmat, role="barrier")] +
                      # +10 * [Layer(width=10e-9, material=Bmat, role="barrier"),
                             # Layer(width=2e-9, material=i_GaAs, role="well"),
                             # Layer(width=7e-9, material=QWmat, role="well"),
                             # Layer(width=2e-9, material=i_GaAs, role="well"),
                             # Layer(width=10e-9, material=Bmat, role="barrier"),] +
                      [Layer(width=2e-9, material=Bmat, role="barrier")] +
                      [Layer(width=50e-9, material=p_GaAs, role="Base")],
                      substrate=i_GaAs)
    print(struc)
    output_2 = QM.schrodinger(struc, quasiconfined=0.05, graphtype='potentialsLDOS', num_eigenvalues=200, show=True)
def showschrodinger2():
    T = 293
    n_GaAs = material("GaAs")(T=T, Nd=1e24)
    p_GaAs = material("GaAs")(T=T, Na=8e22)
    QWmat = material("InGaAs")(T=T, In=0.2, strained=True)
    Bmat = material("GaAsP")(T=T, P=0.1, strained=True)
    i_GaAs = material("GaAs")(T=T)
    QW = PDD.QWunit(
        [
            Layer(width=10e-9, material=Bmat, role="barrier"),
            Layer(width=2e-9, material=i_GaAs, role="well"),
            Layer(width=7e-9, material=QWmat, role="well"),
            Layer(width=2e-9, material=i_GaAs, role="well"),
            Layer(width=10e-9, material=Bmat, role="barrier"),
        ],
        T=T,
        repeat=30,  # สร้างมา30ชั้นเอาไว้กันฝน
        substrate=i_GaAs,  # ใช้intrensic GaAs เป็นฐาน
    )
    QW_list = QW.GetEffectiveQW(wavelengths=wl)
    # for i in QW_list:
    #     print(i)
    QWstr = Structure(
        [Layer(width=50e-9, material=n_GaAs, role='Top_layer')] +
        QW_list +
        [Layer(width=50e-9, material=p_GaAs, role='bottom')], substrate=i_GaAs)
    print(QWstr)
    output_2 = QM.schrodinger(QWstr, quasiconfined=0, graphtype='potentialsLDOS', num_eigenvalues=100, show=True)
def solar_cells():
    T = 298
    wl = np.linspace(350, 1050, 301) * 1e-9

    # สร้างวัสดุทำรับQW ในsolar cell
    QWmat = material("InGaAs")(T=T, In=0.2, strained=True)
    Bmat = material("GaAsP")(T=T, P=0.1, strained=True)
    i_GaAs = material("GaAs")(T=T)
    i_GaAs_dope = material("GaAs")(T=T, Na=1e14)

    # กำหนดโครงสร้างของwellให้มี3ชั้น 1ใช่GaAsP เป็น barrier 2ใช้internsic GaAs เป็นinterlayers 3 ใช้ InGaAs เป็นwell
    QW = PDD.QWunit(
        [
            Layer(width=5e-9, material=Bmat, role="barrier"),
            Layer(width=2e-9, material=i_GaAs_dope, role="well"),
            Layer(width=6e-9, material=QWmat, role="well"),
            Layer(width=2e-9, material=i_GaAs_dope, role="well"),
            Layer(width=5e-9, material=Bmat, role="barrier")
        ],
        T=T,
        repeat=30,  #
        substrate=i_GaAs,
    )

    QW_list = QW.GetEffectiveQW(
        wavelengths=wl)  # optimic ให้QW ดีขึ้น(ไม่แน่ใจทำยังไง)เท่าที่เข้าใจคือมันปรับstrain ในband และปรับค่า absorption
    n_GaAs = material("GaAs")(T=T, Nd=1e24)
    p_GaAs = material("GaAs")(T=T, Na=8e22)

    GaAs_junction = Junction(
    [Layer(width=150e-9, material=n_GaAs, role="Emitter"),]
        + QW_list
        + [Layer(width=2000e-9, material=p_GaAs, role="Base"),],
        sn=1e6,sp=1e6,T=T,kind="PDD",)
    my_solar_cell = SolarCell([GaAs_junction,],
        T=T,substrate=n_GaAs,)
    light_source = LightSource(
        source_type="standard",
        version="AM1.5g",
        x=wl,
        output_units="photon_flux_per_m",
        concentration=1,
    )
    solar_cell_solver(my_solar_cell, "qe",
                      user_options={"light_source": light_source,
                                    "wavelength": wl,
                                    "optics_method": "TMM", }, )

    num_con = 5  # จำนวนในการสร้างแสง
    con_ligth = np.logspace(0, 3, num_con)
    vint = np.linspace(-3.5, 4, 600)
    V = np.linspace(-3.5, 0, 300)
    isc = [];voc = [];FF = [];pmpp = []
    for i, c in enumerate(con_ligth):  # ทำการยิงแสงทั้งหมด20ครั้งตามcon
        light_source.concentration = c
        solar_cell_solver(my_solar_cell, "iv"
                          , user_options={"light_source": light_source,
                                          "wavelength": wl,
                                          "optics_method": None,
                                          "light_iv": True,
                                          "mpp": True,
                                          "voltages": V,
                                          "internal_voltages": vint,
                                          }, )
        isc.append(my_solar_cell.iv["Isc"])
        voc.append(my_solar_cell.iv["Voc"])
        FF.append(my_solar_cell.iv["FF"])
        pmpp.append(my_solar_cell.iv["Pmpp"])
        print('======================================================================')
        print('======================================================================')
        print(f'{(i + num_con ) / (num_con) * 100}%')
        print('======================================================================')
        print('======================================================================')


    with open('data_QW.npy', 'wb') as fout:
        np.save(fout, np.array(isc))
        np.save(fout, np.array(voc))
        np.save(fout, np.array(FF))
        np.save(fout, np.array(pmpp))
        np.save(fout, np.array(con_ligth))
    with open('data_solar_cell.pickle','wb') as out:
        pickle.dump(my_solar_cell, out)
    print('save succeed')
    print('ploting')


    fig2, axes = plt.subplots(2, 2, figsize=(11.25, 8))

    axes[0, 0].semilogx(con_ligth, np.array(pmpp) / con_ligth / 10,"r-o")
    axes[0, 0].set_xlabel("Concentration (suns)")
    axes[0, 0].set_ylabel("Efficiency (%)")

    axes[0, 1].loglog(con_ligth, abs(np.array(isc)), "b-o")
    axes[0, 1].set_xlabel("Concentration (suns)")
    axes[0, 1].set_ylabel("I$_{SC}$ (Am$^{-2}$)")

    axes[1, 0].semilogx(con_ligth, abs(np.array(voc)), "g-o")
    axes[1, 0].set_xlabel("Concentration (suns)")
    axes[1, 0].set_ylabel("V$_{OC}$ (V)")

    axes[1, 1].semilogx(con_ligth, abs(np.array(FF)) * 100, "k-o")
    axes[1, 1].set_xlabel("Concentration (suns)")
    axes[1, 1].set_ylabel("Fill Factor (%)")

    plt.tight_layout()
    fig, ax1 = plt.subplots(1, 1, figsize=(6, 4))

    ax1.plot(wl * 1e9, my_solar_cell.absorbed, "k", label="Total Absorbed")
    ax1.legend(loc="upper right", frameon=False)
    ax1.set_xlabel("Wavelength (nm)")
    ax1.set_ylabel("EQE")
    ax1.set_ylim(0, 1.1)
    ax1.set_xlim(350, 1150)

    plt.tight_layout()
    print('plot finite')
    print('save file')


    fig2.savefig(f'performance_QW_ref.svg', dpi=300)
    fig.savefig(f'EQE_QW_ref.svg', dpi=300)
    text = ''
    data = save_tuple2text(text, 'profile'  + 'width', my_solar_cell)
    save_file_direction(f'final_version_QD', data, f'QD_data')

    def movefile(file, direction):
        save_path = os.path.join(current_path, direction)
        fig1_loc = os.path.join(current_path, file)
        fig1_loc_new = os.path.join(save_path, file)
        shutil.move(fig1_loc, fig1_loc_new)

    current_path = os.getcwd()
    movefile(f'EQE_QW_ref.svg', f'final_version_QD')
    movefile(f'performance_QW_ref.svg', f'final_version_QD')
    # movefile(f'IV_curve_QW_{dot:.2e}_width_{version}.png', f'final_version_QD')
    print('save photo')

def plot_data_from_file(file, fileplk):
    with open(file, 'rb') as filein:
        isc = np.load(filein)
        voc = np.load(filein)
        FF = np.load(filein)
        pmpp = np.load(filein)
        con_ligth = np.load(filein)
    print(isc)
    print(voc)
    print(FF)
    print(pmpp)
    print(con_ligth)

    fig2, axes = plt.subplots(2, 2, figsize=(11.25, 8))

    axes[0, 0].semilogx(con_ligth, np.array(pmpp) / con_ligth / 10, "r-o")
    axes[0, 0].set_xlabel("Concentration (suns)")
    axes[0, 0].set_ylabel("Efficiency (%)")

    axes[0, 1].loglog(con_ligth, abs(np.array(isc)/10000), "b-o")
    axes[0, 1].set_xlabel("Concentration (suns)")
    axes[0, 1].set_ylabel("I$_{SC}$ (Am$^{-2}$)")

    axes[1, 0].semilogx(con_ligth, abs(np.array(voc)), "g-o")
    axes[1, 0].set_xlabel("Concentration (suns)")
    axes[1, 0].set_ylabel("V$_{OC}$ (V)")

    axes[1, 1].semilogx(con_ligth, abs(np.array(FF)) * 100, "k-o")
    axes[1, 1].set_xlabel("Concentration (suns)")
    axes[1, 1].set_ylabel("Fill Factor (%)")

    plt.tight_layout()
    # fig, ax1 = plt.subplots(1, 1, figsize=(6, 4))
    # ax1.plot(wl * 1e9, my_solar_cell.absorbed, "k", label="Total Absorbed")
    # ax1.legend(loc="upper right", frameon=False)
    # ax1.set_xlabel("Wavelength (nm)")
    # ax1.set_ylabel("EQE")
    # ax1.set_ylim(0, 1.1)
    # ax1.set_xlim(350, 1150)
    # plt.tight_layout()
    print('plot finite')
    print('save file')
    fig2.savefig(f'performance_QW_ref.svg', dpi=300)
    plt.show()
    # fig.savefig(f'EQE_QW_ref.svg', dpi=300)
    # text = ''
    # data = save_tuple2text(text, 'profile'  + 'width', my_solar_cell)
    # save_file_direction(f'final_version_QD', data, f'QD_data')

    def movefile(file, direction):
        save_path = os.path.join(current_path, direction)
        fig1_loc = os.path.join(current_path, file)
        fig1_loc_new = os.path.join(save_path, file)
        shutil.move(fig1_loc, fig1_loc_new)

    current_path = os.getcwd()
    # movefile(f'EQE_QW_ref.svg', f'final_version_QD')
    movefile(f'performance_QW_ref.svg', f'final_version_QD')
    # movefile(f'IV_curve_QW_{dot:.2e}_width_{version}.png', f'final_version_QD')
    print(current_path, 'path direction')
    print('save photo')
    # print(my_solar_cell)
    print('succeed')
# solar_cells()
plot_data_from_file("data_QW.npy","data_solar_cell.pickle")