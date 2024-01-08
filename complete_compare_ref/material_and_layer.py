from solcore import si, material
from solcore.structure import Layer, Structure, Junction
from solcore.solar_cell import SolarCell
import numpy as np

#===============================================================
#setup
T = 300
def printstructure(solar_cell):
    space = "+" +'='*70 + "+"
    print(space+ '\n' + f"\n".join(["{}".format(layer) for layer in solar_cell])+ '\n' + space)


#===============================================================
#===============================================================
#material (ref : High-efficiency GaAs and GaInP solar cells grown by all solid-state Molecular-Beam-Epitaxy)

n_GaAs_window = material("GaAs")(T=T, Nd=si('4e18 cm-3'))
n_AlInP = material("AlInP")(T=T, Al=0.42, Nd=si('3e17 cm-3'))
n_GaAs = material("GaAs")(T=T, Nd=si("2e18 cm-3"))
p_GaAs = material("GaAs")(T=T, Na=si("5e16 cm-3"))
p_GaInP = material("GaInP")(T=T, In=0.35, Na=si("2e18 cm-3"))
p_GaAs_buffer = material("GaAs")(T=T, Na=si("2e18 cm-3"))
p_GaAs_sub = material('GaAs')(T=T, Na=si('2e17 cm-3'))
#layer

GaAs_junction = Junction([
    Layer(width=si("300 nm"), material=n_GaAs_window, role="front_contact"),
    Layer(width=si("30 nm"), material=n_AlInP, role="window"),
    Layer(width=si("100 nm"), material=n_GaAs, role="Emitter"),
    Layer(width=si("2800 nm"), material=p_GaAs, role="Base"),
    Layer(width=si("100 nm"), material=p_GaInP, role="BSF"),
    Layer(width=si("150 nm"), material=p_GaAs_buffer, role="Buffer"),

    ],
    T=T, kind="PDD", )
my_solar_cell = SolarCell([

    GaAs_junction,

    ]
    ,T=T,substrate=p_GaAs_sub,)
# printstructure(my_solar_cell)

#===============================================================
#vary GaInP
vary_GaInP = list()
con_GaInP = np.linspace(0.2,0.6,10)
for i in con_GaInP:
    n_GaAs_window = material("GaAs")(T=T, Nd=si('4e18 cm-3'))
    n_AlInP = material("AlInP")(T=T, Al=0.32, Nd=si('3e17 cm-3'))
    n_GaAs = material("GaAs")(T=T, Nd=si("2e18 cm-3"))
    p_GaAs = material("GaAs")(T=T, Na=si("5e16 cm-3"))
    p_GaInP = material("GaInP")(T=T, In=i, Na=si("2e18 cm-3"))
    p_GaAs_buffer = material("GaAs")(T=T, Na=si("2e18 cm-3"))
    p_GaAs_sub = material('GaAs')(T=T, Na=si('2e17 cm-3'))
    #layer

    GaAs_junction = Junction([
        Layer(width=si("300 nm"), material=n_GaAs_window, role="front_contact"),
        Layer(width=si("30 nm"), material=n_AlInP, role="window"),
        Layer(width=si("100 nm"), material=n_GaAs, role="Emitter"),
        Layer(width=si("2800 nm"), material=p_GaAs, role="Base"),
        Layer(width=si("100 nm"), material=p_GaInP, role="BSF"),
        Layer(width=si("150 nm"), material=p_GaAs_buffer, role="Buffer"),
        ],
        T=T, kind="PDD", )
    vary_GaInP.append(SolarCell([

        GaAs_junction,

        ]
        ,T=T,substrate=p_GaAs_sub,))
#===============================================================
#vary AlInP
vary_AlInP = list()
con_AlInP = np.linspace(0.1,0.7,10)
for i in con_AlInP:
    n_GaAs_window = material("GaAs")(T=T, Nd=si('4e18 cm-3'))
    n_AlInP = material("AlInP")(T=T, Al=i, Nd=si('3e17 cm-3'))
    n_GaAs = material("GaAs")(T=T, Nd=si("2e18 cm-3"))
    p_GaAs = material("GaAs")(T=T, Na=si("5e16 cm-3"))
    p_GaInP = material("GaInP")(T=T, In=0.2, Na=si("2e18 cm-3"))
    p_GaAs_buffer = material("GaAs")(T=T, Na=si("2e18 cm-3"))
    p_GaAs_sub = material('GaAs')(T=T, Na=si('2e17 cm-3'))
    #layer

    GaAs_junction = Junction([
        Layer(width=si("300 nm"), material=n_GaAs_window, role="front_contact"),
        Layer(width=si("30 nm"), material=n_AlInP, role="window"),
        Layer(width=si("100 nm"), material=n_GaAs, role="Emitter"),
        Layer(width=si("2800 nm"), material=p_GaAs, role="Base"),
        Layer(width=si("100 nm"), material=p_GaInP, role="BSF"),
        Layer(width=si("150 nm"), material=p_GaAs_buffer, role="Buffer"),
        ],
        T=T, kind="PDD", )
    vary_AlInP.append(SolarCell([

        GaAs_junction,

        ]
        ,T=T,substrate=p_GaAs_sub,))
#===============================================================
vary_AlInP_and_GaInP = []
buffer = []
for i in con_AlInP:
    for j in con_GaInP:
        n_GaAs_window = material("GaAs")(T=T, Nd=si('4e18 cm-3'))
        n_AlInP = material("AlInP")(T=T, Al=i, Nd=si('3e17 cm-3'))
        n_GaAs = material("GaAs")(T=T, Nd=si("2e18 cm-3"))
        p_GaAs = material("GaAs")(T=T, Na=si("5e16 cm-3"))
        p_GaInP = material("GaInP")(T=T, In=j, Na=si("2e18 cm-3"))
        p_GaAs_buffer = material("GaAs")(T=T, Na=si("2e18 cm-3"))
        p_GaAs_sub = material('GaAs')(T=T, Na=si('2e17 cm-3'))
        #layer

        GaAs_junction = Junction([
            Layer(width=si("300 nm"), material=n_GaAs_window, role="front_contact"),
            Layer(width=si("30 nm"), material=n_AlInP, role="window"),
            Layer(width=si("100 nm"), material=n_GaAs, role="Emitter"),
            Layer(width=si("2800 nm"), material=p_GaAs, role="Base"),
            Layer(width=si("100 nm"), material=p_GaInP, role="BSF"),
            Layer(width=si("150 nm"), material=p_GaAs_buffer, role="Buffer"),
            ],
            T=T, kind="PDD", )
        buffer.append(SolarCell([

            GaAs_junction,

            ]
            ,T=T,substrate=p_GaAs_sub,))
    vary_AlInP_and_GaInP.append(buffer)
    buffer = []
#===============================================================
#arcoding
solat_cell_with_arc = []
n_GaAs_window = material("GaAs")(T=T, Nd=si('4e18 cm-3'))
n_AlInP = material("AlInP")(T=T, Al=0.42, Nd=si('3e17 cm-3'))
n_GaAs = material("GaAs")(T=T, Nd=si("2e18 cm-3"))
p_GaAs = material("GaAs")(T=T, Na=si("5e16 cm-3"))
p_GaInP = material("GaInP")(T=T, In=0.35, Na=si("2e18 cm-3"))
p_GaAs_buffer = material("GaAs")(T=T, Na=si("2e18 cm-3"))
p_GaAs_sub = material('GaAs')(T=T, Na=si('2e17 cm-3'))
#layer

GaAs_junction = Junction([
    Layer(width=si("300 nm"), material=n_GaAs_window, role="front_contact"),
    Layer(width=si("30 nm"), material=n_AlInP, role="window"),
    Layer(width=si("100 nm"), material=n_GaAs, role="Emitter"),
    Layer(width=si("2800 nm"), material=p_GaAs, role="Base"),
    Layer(width=si("100 nm"), material=p_GaInP, role="BSF"),
    Layer(width=si("150 nm"), material=p_GaAs_buffer, role="Buffer"),

    ],
    T=T, kind="PDD", )

SiO2 = material("SiO2")()
Si3N4 = material("Si3N4")()
MgF2 = material("MgF2")()
ZnS = material("ZnScub")()
solat_cell_with_arc = SolarCell([
    # Layer(width=si("402.87 nm"), material=SiO2, role="ARC1"),
    # Layer(width=si("287 nm"), material=Si3N4, role="ARC2"),
    Layer(width=si("110 nm"), material=MgF2, role="ARC1"),
    Layer(width=si("60 nm"), material=ZnS , role="ARC2"),

    GaAs_junction,

    ]
    ,T=T,substrate=p_GaAs_sub,)
printstructure(solat_cell_with_arc)