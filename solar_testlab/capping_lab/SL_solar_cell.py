from solcore import si, material
from solcore.structure import Layer, Structure, Junction, TunnelJunction
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
wl = np.linspace(350, 1200, 401) * 1e-9

light_source = LightSource(
    source_type="standard",
    version="AM1.5g",
    x=wl,
    output_units="photon_flux_per_m",
    concentration=1,
)


def solar_cells(version,active_dot=True):

    T = 298
    wl = np.linspace(350, 1200, 301) * 1e-9

    # First, we create the materials of the QW
    # สร้างวัสดุทำรับQW ในsolar cell
    n_GaAs = material('GaAs')(T=T, Nd=1e22)
    p_GaAs = material('GaAs')(T=T, Na=8e20)
    barrier = material('AlGaAs')(T=T, Al=0.3, Nd=1e22)
    QWmat1 = material("InSb")(T=T, strained=True, electron_mobility = 7.7 , hole_mobility=0.0850)
    QWmat2 = material("GaSb")(T=T, strained=True, electron_mobility= 0.3, hole_mobility= 0.1)
    i_GaAs = material("GaAs")(T=T)
    if active_dot:
        QW = PDD.QWunit([
        Layer(width=100e-9, material=barrier, role='barrier'),
        Layer(width=83.84e-9, material=i_GaAs, role="interlayer"),
        Layer(width=16.16e-9, material=QWmat1, role="well"),
        Layer(width=96.46e-9, material=i_GaAs, role="interlayer"),
        Layer(width=2.54e-9, material=QWmat2, role="well"),
        Layer(width=50e-9, material=n_GaAs, role="interlayer"),
        Layer(width=100e-9, material=barrier, role="barrier"),
                        ],substrate=p_GaAs, T=T, repeat=1)
        QW_list = QW.GetEffectiveQW(wavelengths=wl)
        #

        GaAs_junction = Junction(
            [
                Layer(width=40e-9, material=n_GaAs, role="window"),
                Layer(width=10e-9, material=n_GaAs, role="window"),
            ]
             # Layer(width=barrier, material=Bmat, role="barrier")]
            + QW_list
            # + [Layer(width=barrier, material=Bmat, role="barrier"),
            +[
                Layer(width=100e-9, material=n_GaAs, role="Emitter"),
                Layer(width=200e-9, material=p_GaAs, role='buffer'),
                # Layer(width=2000e-9, material=p_GaAs, role="Base"),
            ],T=T,kind="PDD",)
        my_solar_cell = SolarCell(
            [
            GaAs_junction,
            ]
            , T=T, substrate=p_GaAs, )
    else:
        n_GaAs = material('GaAs')(T=T, Nd=1e22)
        p_GaAs = material('GaAs')(T=T, Na=8e20)
        barrier = material('AlGaAs')(T=T, Al=0.3, Nd=1e20)
        QWmat1 = material("InSb")(T=T, strained=True, electron_mobility=7.7, hole_mobility=0.0850)
        QWmat2 = material("GaSb")(T=T, strained=True, electron_mobility=0.3, hole_mobility=0.1)
        i_GaAs = material("GaAs")(T=T)
        GaAs_junction = Junction(
            [
        Layer(width=100e-9, material=n_GaAs, role="Emitter"),
        Layer(width=200e-9, material=i_GaAs, role='buffer'),
            # Layer(width=barrier, material=Bmat, role="barrier")]
            # + QW_list
                # + [Layer(width=barrier, material=Bmat, role="barrier"),
        Layer(width=2000e-9, material=p_GaAs, role="Base"),
            ],
            sn=1e6,
            sp=1e6, T=T, kind="PDD", )
        my_solar_cell = SolarCell([Layer(width=50e-9, material=n_GaAs, role="Emitter"),
        Layer(width=100e-9, material=barrier, role='barrier'),
        Layer(width=250-9, material=n_GaAs, role="interlayer"),
        Layer(width=100e-9, material=barrier, role="barrier"),
            GaAs_junction, ]
            , T=T, substrate=p_GaAs, )


    solar_cell_solver(my_solar_cell, "qe",
                      user_options={"light_source": light_source,
                                    "wavelength": wl,
                                    "optics_method": "TMM",}, )

    num_con = 7 # จำนวนในการสร้างแสง
    con = np.logspace(0, 3, num_con)
    vint = np.linspace(-3.5, 4, 600)
    V = np.linspace(-3.5, 0, 300)
    allI = []; isc = []; voc = []; FF = []; pmpp = []
    fig3, axIV = plt.subplots(1, 1, figsize=(6, 4))
    fig4, axJV = plt.subplots(1, 1, figsize=(6, 4))
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
        axJV.semilogy(my_solar_cell[0].voltage, abs(my_solar_cell[0].current), 'k', linewidth=4, label='Total')
        axJV.semilogy(my_solar_cell[0].voltage, abs(my_solar_cell[0].recombination_currents['Jrad']), 'r', label='Jrad')
        axJV.semilogy(my_solar_cell[0].voltage, abs(my_solar_cell[0].recombination_currents['Jsrh']), 'b', label='Jsrh')
        axJV.semilogy(my_solar_cell[0].voltage, abs(my_solar_cell[0].recombination_currents['Jsur']), 'g', label='Jsur')
        if active_dot:
            print('======================================================================')
            print('======================================================================')
            print(f'{(i)/(num_con)*100}%')
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
    axIV.set_xlim(0, 2)
    axIV.set_xlabel("Voltage (V)")
    axIV.set_ylabel("Normalized current (-)")
    plt.tight_layout()
    axJV.legend()
    axJV.set_xlim(-0.5, 1.3)
    axJV.set_ylim(1e-10, 1e5)
    axJV.set_xlabel('Bias (V)')
    axJV.set_ylabel('Current (A/m$^2}$)')
    plt.tight_layout()


    fig2, axes = plt.subplots(2, 2, figsize=(11.25, 8))

    axes[0, 0].semilogx(con, np.array(pmpp) / con / 10,"r-o")
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
    fig0, ax1 = plt.subplots(1, 1, figsize=(6, 4))
    for j in my_solar_cell.junction_indices:  # junctionคือหยั่ง
        zz = (my_solar_cell[j].short_circuit_data.Bandstructure["x"] + my_solar_cell[j].offset)
        n = my_solar_cell[j].short_circuit_data.Bandstructure["n"]
        p = my_solar_cell[j].short_circuit_data.Bandstructure["p"]
        ax1.semilogy(zz * 1e9, n, "b")  # อันนี้น่าจะเป็นการระบุcarrier densityของตอนฉายแสงของe และhole
        ax1.semilogy(zz * 1e9, p, "r")

        zz = my_solar_cell[j].equilibrium_data.Bandstructure["x"] + my_solar_cell[j].offset
        n = my_solar_cell[j].equilibrium_data.Bandstructure["n"]
        p = my_solar_cell[j].equilibrium_data.Bandstructure["p"]
        ax1.semilogy(zz * 1e9, n, "b--")  # อันนี้น่าจะเป็นการระบุcarrier densityของตอนปิดแสงของe และhole
        ax1.semilogy(zz * 1e9, p, "r--")
    plt.tight_layout()

    ax1.set_xlabel("Position (nm)")
    ax1.set_ylabel("Carrier density (m$^{-3}$)")
    plt.tight_layout()
    fig3.savefig(f'IV_curve_QW_{version}.png', dpi=300)
    fig2.savefig(f'performance_QW_{version}.png', dpi=300)
    fig.savefig(f'EQE_QW_{version}.png', dpi=300)
    fig0.savefig(f'carrier_distribution{version}.png', dpi=300)
    text = ''
    data = save_tuple2text(text,'profile'+version, my_solar_cell)
    for x in my_solar_cell:
        data += str(x)
    save_file_direction(f'data_of_{version}', data, f'QW{version}')
    def movefile(file,direction):
        save_path = os.path.join(current_path, direction)
        fig1_loc = os.path.join(current_path, file)
        fig1_loc_new = os.path.join(save_path, file)
        shutil.move(fig1_loc, fig1_loc_new)

    current_path = os.getcwd()
    movefile(f'EQE_QW_{version}.png', f'data_of_{version}')
    movefile(f'performance_QW_{version}.png', f'data_of_{version}')
    movefile(f'IV_curve_QW_{version}.png', f'data_of_{version}')
    movefile(f'carrier_distribution{version}.png', f'data_of_{version}')

    plt.show()
solar_cells('SL_version_3_equaldopind_barrier', active_dot=True)