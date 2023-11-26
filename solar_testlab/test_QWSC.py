from solcore import material
from solcore.structure import Layer, Junction, TunnelJunction
from solcore.solar_cell import SolarCell
from solcore.solar_cell_solver import solar_cell_solver
from solcore.light_source import LightSource
import solcore.poisson_drift_diffusion as PDD

import numpy as np
import matplotlib.pyplot as plt
import pickle
list_ice = [int(i*10)for i in range(1,6)]
for layers in list_ice:
    T = 298
    wl = np.linspace(350, 1050, 301) * 1e-9

    # First, we create the materials of the QW
    # สร้างวัสดุทำรับQW ในsolar cell
    QWmat = material("InGaAs")(T=T, In=0.2, strained=True)
    Bmat = material("GaAsP")(T=T, P=0.1, strained=True)
    i_GaAs = material("GaAs")(T=T)
    i_GaAs_dope = material("GaAs")(T=T, Na= 1e14)

    # The QW is 7 nm wide, with GaAs interlayers 2 nm thick at each side and GaAsP barriers
    # 10 nm thick. The final device will have 30 of these QWs.
    # กำหนดโครงสร้างของwellให้มี3ชั้น 1ใช่GaAsP เป็น barrier 2ใช้internsic GaAs เป็นinterlayers 3 ใช้ InGaAs เป็นwell
    QW = PDD.QWunit(
        [
             Layer(width=10e-9, material=Bmat, role="barrier"),
             Layer(width=2e-9, material=i_GaAs, role="well"),
             Layer(width=ุ6e-9, material=QWmat, role="well"),
             Layer(width=2e-9, material=i_GaAs, role="well"),
             Layer(width=10e-9, material=Bmat, role="barrier")
        ],
        T=T,
        repeat=layers, #สร้างมา30ชั้นเอาไว้กันฝน
        substrate=i_GaAs, #ใช้intrensic GaAs เป็นฐาน
    )

    # We solve the quantum properties of the QW, leaving the default values of all
    # parameters
    QW_list = QW.GetEffectiveQW(wavelengths=wl) #optimic ให้QW ดีขึ้น(ไม่แน่ใจทำยังไง)เท่าที่เข้าใจคือมันปรับstrain ในband และปรับค่า absorption
    #แบบ่ว่าเบิ่้มแต่ไม่รู้หลักการเท่าไร

    # Materials for the BOTTOM junction
    # ทำส่วนล่าง
    window_bottom = material("GaInP")(T=T, Nd=5e24, In=0.49)
    n_GaAs = material("GaAs")(T=T, Nd=1e24)
    p_GaAs = material("GaAs")(T=T, Na=8e22)
    bsf_bottom = material("GaInP")(T=T, Na=5e24, In=0.49)

    # If you want to test the code without QWs, to make ti a bit faster, comment the line
    # with QW_list ประกอบร่างให้กับsolar cell
    GaAs_junction = Junction(
        [
            Layer(width=10e-9, material=window_bottom, role="Window"),#ชั้นบน
            Layer(width=150e-9, material=n_GaAs, role="Emitter"),
            # Layer(width=2e-9, material=Bmat, role="barrier")
    ]
        # Comment the following line to remove the QWs
        # + QW_list
        + [
            # Layer(width=2e-9, material=Bmat, role="barrier"),
            Layer(width=2000e-9, material=p_GaAs, role="Base"),
            Layer(width=200e-9, material=bsf_bottom, role="BSF"),#ชั้นล่าง

        ],
        sn=1e6,
        sp=1e6,
        T=T,
        kind="PDD",
    )
    # And the materials needed for the anti reflecting coating
    MgF2 = material("MgF2")()
    ZnS = material("ZnScub")()
    # Finally, we put everithing together to make a solar cell
    my_solar_cell = SolarCell(
        [
            Layer(width=110e-9, material=MgF2, role="ARC1"),
            Layer(width=60e-9, material=ZnS, role="ARC2"),
            GaAs_junction,
        ],
        T=T,
        substrate=n_GaAs,
    )
    light_source = LightSource(
        source_type="standard",
        version="AM1.5g",
        x=wl,
        output_units="photon_flux_per_m",
        concentration=1,
    )

    # The definitions are all done, so we just start solving the properties,
    # starting with the QE. We calculate the QE curve under illumination
    # สร้างเสร็จแล้วก็แก้สมการได้เย้ๆ งงอะไรคือQE(ใช่่ quantum efficiency หรือป่าว)
    solar_cell_solver(my_solar_cell,"qe",
        user_options={"light_source": light_source,
            "wavelength": wl,
            "optics_method": "TMM",
        },)

    # And now, the IV curves under various concentration levels.
    # NOTE: Due to the presence of QWs and the fact we calculate things a 19 different
    # concentrations, this might take a while (~4 hours).
    # Remove the QWs as indicated above to test the code much faster.

    num_con = 5 #จำนวนในการสร้างแสง
    con = np.logspace(0, 3, num_con)
    vint = np.linspace(-3.5, 4, 600)
    V = np.linspace(-3.5, 0, 300)

    allI = []
    isc = []
    voc = []
    FF = []
    pmpp = []

    fig3, axIV = plt.subplots(1, 1, figsize=(6, 4))
    for c in con:#ทำการยิงแสงทั้งหมด20ครั้งตามcon
        light_source.concentration = c
        solar_cell_solver(
            my_solar_cell,#วัสดุ
            "iv",#วัดIV
            user_options={
                "light_source": light_source,
                "wavelength": wl,
                "optics_method": None,#ไม่คำนึงถึงการตกกระทบทางแสงหรือเข้า100%ไม่แน่ใจ
                "light_iv": True,#ให่IV ออกมา
                "mpp": True,#ให้Mppออกมา
                "voltages": V,#ให้Vออกมา
                "internal_voltages": vint,#ไม่แน่ใจอะไรคือVintเป็นbulid in หรือป่าวเท้าที่เห็นเป็นVoc
            },
        )

        isc.append(my_solar_cell.iv["Isc"]) #นำค่าIscที่ได้มาในแต่ละความเข้มแสงมาเก็บไว้
        voc.append(my_solar_cell.iv["Voc"]) #นำค่าVocที่ได้มานั้นมาเก็บ
        FF.append(my_solar_cell.iv["FF"])   #FF(Isc*Vsc/Pmpp)
        pmpp.append(my_solar_cell.iv["Pmpp"])#Pmpp power maximum power peak
        allI.append(my_solar_cell.iv["IV"][1])#ไม่แ่นใจ่ว่าทำไมต้อง1แต่คิดว่า0น่าจะเป็นแรงดันต้องลองรันในspyder(ลองดูสมบัติของmy_solar_vell.iv['IV][1]

        # outpath = "C:\\Users\\usuwa\\Desktop\\senior_project\\simulation_solar_cell"
        # fname = f'light_conc_{c}.pkl'
        # with open(outpath+'\\'+fname,'wb') as fout:
        #     pickle.save(my_solar_cell,fout)

        # And now, everything is plotting...
        axIV.plot(-V, my_solar_cell.iv["IV"][1] / isc[-1], label=int(c))#เป็นการวัดลักษณะIVเมิอแสงเปลี่ยนไป

    axIV.legend(loc="lower left", frameon=False)
    axIV.set_ylim(0, 1.1)
    axIV.set_xlim(0, 3.5)
    axIV.set_xlabel("Voltage (V)")
    axIV.set_ylabel("Normalized current (-)")

    plt.tight_layout()
    fig3.savefig(f'IV_curve_QW_{layers}_layers.svg', dpi=300)

    fig2, axes = plt.subplots(2, 2, figsize=(11.25, 8))

    axes[0, 0].semilogx(con, np.array(pmpp) / con / 10, "r-o") #ทำไมต้องหาร10ด้วยไม่เข้าใจเลย(ไม่ใช้หารด้วย 0.97หรอเพราะว่าAM1.5g = 970w/m^2)แล้วทำไมไม่x100
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
    fig2.savefig(f'performance_QW_{layers}_layers.svg', dpi=300)

    # We can plot the electron and hole densities in equilibrium and at short circuit
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.25, 4))
    for j in my_solar_cell.junction_indices:#junctionคือหยั่ง
        zz = (
            my_solar_cell[j].short_circuit_data.Bandstructure["x"] + my_solar_cell[j].offset)
        n = my_solar_cell[j].short_circuit_data.Bandstructure["n"]
        p = my_solar_cell[j].short_circuit_data.Bandstructure["p"]
        ax1.semilogy(zz * 1e9, n, "b")#อันนี้น่าจะเป็นการระบุcarrier densityของตอนฉายแสงของe และhole
        ax1.semilogy(zz * 1e9, p, "r")

        zz = my_solar_cell[j].equilibrium_data.Bandstructure["x"] + my_solar_cell[j].offset
        n = my_solar_cell[j].equilibrium_data.Bandstructure["n"]
        p = my_solar_cell[j].equilibrium_data.Bandstructure["p"]
        ax1.semilogy(zz * 1e9, n, "b--")#อันนี้น่าจะเป็นการระบุcarrier densityของตอนปิดแสงของe และhole
        ax1.semilogy(zz * 1e9, p, "r--")

    ax1.set_xlabel("Position (nm)")
    ax1.set_ylabel("Carrier density (m$^{-3}$)")
    plt.tight_layout()

    # And we plot the QE
    labels = ["EQE GaInP", "EQE GaAs"]
    colours = ["b", "r"]
    for i, j in enumerate(my_solar_cell.junction_indices):#เป็นการพล็อตของeqeของแต่ละอันของGaAs และGaInP
        ax2.plot(wl * 1e9, my_solar_cell[j].eqe(wl), colours[i], label=labels[i])

    ax2.plot(wl * 1e9, my_solar_cell.absorbed, "k", label="Total Absorbed")

    ax2.legend(loc="upper right", frameon=False)
    ax2.set_xlabel("Wavelength (nm)")
    ax2.set_ylabel("EQE")
    ax2.set_ylim(0, 1.1)
    ax2.set_xlim(350, 1150)

    plt.tight_layout()
    fig.savefig(f'EQE_QW_{layers}_layers.svg', dpi=300)

# plt.show()
