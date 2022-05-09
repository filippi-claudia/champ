    !> This module provides the periodic table data for elements upto Z=18
    !! @author Ravindra Shinde
    !! @date June 23 2021
    !! @remarks The elemental data is taken from the NIST website
    !! (https://physics.nist.gov/PhysRefData/Handbook/periodictable_a.htm)
    !! The data tabulated in this file is only for the most abundant isotope only.

module periodic_table
    implicit none
    public :: element, atom_t
    private :: i1, sp

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: i1 = kind(1)

    type :: atom_t
        real(kind=sp)        :: atomic_mass  ! in amu
        character(len=20)    :: name
        character(len=3)     :: symbol
        integer(kind=i1)     :: znuclear
        integer(kind=i1)     :: core
        integer(kind=i1)     :: nvalence
        integer(kind=i1)     :: isotope = 1  ! currently supports only the most abundant
!   future plans
!                               covalent_radii
!                               ground state term symbol
!                               total spin
!                               magnetic moment

    end type atom_t

    contains

    type(atom_t) function element(sym) result(atom)

    character(len=*) :: sym

    select case(sym)
        case("ghost", "Ghost Atom", "X", "0")

            atom = atom_t(name="ghost",  symbol="X",     atomic_mass=0.0,   znuclear=0, core=0, nvalence=0)

        case("hydrogen", "Hydrogen", "H", "1")

            atom = atom_t(name="hydrogen",  symbol="H",     atomic_mass=1.007825,   znuclear=1, core=0, nvalence=1)

        case("helium", "Helium", "He", "2")

            atom = atom_t(name="helium",    symbol="He",    atomic_mass=4.00260,    znuclear=2, core=0, nvalence=2)

        case("lithium", "Lithium", "Li", "3")

            atom = atom_t(name="lithium",   symbol="Li",    atomic_mass=7.016003,   znuclear=3, core=2, nvalence=1)

        case("beryllium", "Beryllium", "Be", "4")

            atom = atom_t(name="beryllium", symbol="Be",    atomic_mass=9.012182,   znuclear=4, core=2, nvalence=2)

        case("boron", "Boron", "B", "5")

            atom = atom_t(name="boron",     symbol="B",     atomic_mass=11.009305,  znuclear=5, core=2, nvalence=3)

        case("carbon", "Carbon", "C", "6")

            atom = atom_t(name="carbon",    symbol="C",     atomic_mass=12.000000,  znuclear=6, core=2, nvalence=4)

        case("nitrogen", "Nitrogen", "N", "7")

            atom = atom_t(name="nitrogen",  symbol="N",     atomic_mass=14.003074,  znuclear=7, core=2, nvalence=5)

        case("oxygen", "Oxygen", "O", "8")

            atom = atom_t(name="oxygen",    symbol="O",     atomic_mass=15.994915,  znuclear=8, core=2, nvalence=6)

        case("fluorine", "Fluorine", "F", "9")

            atom = atom_t(name="fluorine",  symbol="F",     atomic_mass=18.9984032, znuclear=9, core=2, nvalence=7)

        case("neon", "Neon", "Ne", "10")

            atom = atom_t(name="neon",      symbol="Ne",    atomic_mass=19.992435,  znuclear=10, core=2, nvalence=8)

        case("sodium", "Sodium", "Na", "11")

            atom = atom_t(name="sodium",    symbol="Na",    atomic_mass=22.989767,  znuclear=11, core=10, nvalence=1)

        case("magnesium", "Magnesium", "Mg", "12")

            atom = atom_t(name="magnesium", symbol="Mg",    atomic_mass=23.985042,  znuclear=12, core=10, nvalence=2)

        case("aluminum", "Aluminum", "aluminium", "Aluminium", "Al", "13")

            atom = atom_t(name="aluminum",  symbol="Al",    atomic_mass=26.981540,  znuclear=13, core=10, nvalence=3)

        case("silicon", "Silicon", "Si", "14")

            atom = atom_t(name="silicon",   symbol="Si",    atomic_mass=27.976927,  znuclear=14, core=10, nvalence=4)

        case("phosphorus", "Phosphorus", "P", "15")

            atom = atom_t(name="phosphorus",symbol="P",     atomic_mass=30.973762,  znuclear=15, core=10, nvalence=5)

        case("sulfur", "Sulfur", "sulphur", "Sulphur", "S", "16")

            atom = atom_t(name="sulfur",    symbol="S",     atomic_mass=31.972070,  znuclear=16, core=10, nvalence=6)

        case("chlorine", "Chlorine", "Cl", "17")

            atom = atom_t(name="chlorine",  symbol="Cl",    atomic_mass=34.968852,  znuclear=17, core=10, nvalence=7)

        case("argon", "Argon", "Ar", "18")

            atom = atom_t(name="argon",     symbol="Ar",    atomic_mass=39.962384,  znuclear=18, core=10, nvalence=8)

        case("potassium", "Potassium", "K", "19")

            atom = atom_t(name="potassium", symbol="K",     atomic_mass=38.963707,  znuclear=19, core=18, nvalence=1)

        case("calcium", "Calcium", "Ca", "20")

            atom = atom_t(name="calcium",   symbol="Ca",    atomic_mass=39.962591,  znuclear=20, core=18, nvalence=2)

        case("scandium", "Scandium", "Sc", "21")

            atom = atom_t(name="scandium",  symbol="Sc",    atomic_mass=44.955910,  znuclear=21, core=18, nvalence=3)

        case("titanium", "Titanium", "Ti", "22")

            atom = atom_t(name="titanium",  symbol="Ti",    atomic_mass=47.947947,  znuclear=22, core=18, nvalence=4)

        case("vanadium", "Vanadium", "V", "23")

            atom = atom_t(name="vanadium",  symbol="V",     atomic_mass=50.943962,  znuclear=23, core=18, nvalence=5)

        case("chromium", "Chromium", "Cr", "24")

            atom = atom_t(name="chromium",  symbol="Cr",    atomic_mass=51.940509,  znuclear=24, core=18, nvalence=6)

        case("manganese", "Manganese", "Mn", "25")

            atom = atom_t(name="manganese", symbol="Mn",    atomic_mass=54.938047,  znuclear=25, core=18, nvalence=7)

        case("iron", "Iron", "Fe", "26")

            atom = atom_t(name="iron",      symbol="Fe",    atomic_mass=55.934939,  znuclear=26, core=18, nvalence=8)

        case("cobalt", "Cobalt", "Co", "27")

            atom = atom_t(name="cobalt",    symbol="Co",    atomic_mass=58.933198,  znuclear=27, core=18, nvalence=9)

        case("nickel", "Nickel", "Ni", "28")

            atom = atom_t(name="nickel",    symbol="Ni",    atomic_mass=57.935346,  znuclear=28, core=18, nvalence=10)

        case("copper", "Copper", "Cu", "29")

            atom = atom_t(name="copper",    symbol="Cu",    atomic_mass=62.939598,  znuclear=29, core=18, nvalence=11)

        case("zinc", "Zinc", "Zn", "30")

            atom = atom_t(name="zinc",      symbol="Zn",    atomic_mass=63.929145,  znuclear=30, core=18, nvalence=12)

        case("gallium", "Gallium", "Ga", "31")

            atom = atom_t(name="gallium",   symbol="Ga",    atomic_mass=68.925580,  znuclear=31, core=18, nvalence=13)

        case("germanium", "Germanium", "Ge", "32")

            atom = atom_t(name="germanium", symbol="Ge",    atomic_mass=73.921177,  znuclear=32, core=18, nvalence=14)

        case("arsenic", "Arsenic", "As", "33")

            atom = atom_t(name="arsenic",   symbol="As",    atomic_mass=74.921594,  znuclear=33, core=18, nvalence=15)

        case("selenium", "Selenium", "Se", "34")

            atom = atom_t(name="selenium",  symbol="Se",    atomic_mass=79.916520,  znuclear=34, core=18, nvalence=16)

        case("bromine", "Bromine", "Br", "35")

            atom = atom_t(name="bromine",   symbol="Br",    atomic_mass=78.918336,  znuclear=35, core=18, nvalence=17)

        case("krypton", "Krypton", "Kr", "36")

            atom = atom_t(name="krypton",   symbol="Kr",    atomic_mass=83.911507,  znuclear=36, core=18, nvalence=18)

        case("rubidium", "Rubidium", "Rb", "37")

            atom = atom_t(name="rubidium",  symbol="Rb",    atomic_mass=84.911794,  znuclear=37, core=36, nvalence=1)

        case("strontium", "Strontium", "Sr", "38")

            atom = atom_t(name="strontium", symbol="Sr",    atomic_mass=87.905619,  znuclear=38, core=36, nvalence=2)

        case("yttrium", "Yttrium", "Y", "39")

            atom = atom_t(name="yttrium",   symbol="Y",     atomic_mass=88.905849,  znuclear=39, core=36, nvalence=3)

        case("zirconium", "Zirconium", "Zr", "40")

            atom = atom_t(name="zirconium", symbol="Zr",    atomic_mass=89.904703,  znuclear=40, core=36, nvalence=4)

        case("niobium", "Niobium", "Nb", "41")

            atom = atom_t(name="niobium",   symbol="Nb",    atomic_mass=92.906377,  znuclear=41, core=36, nvalence=5)

        case("molybdenum", "Molybdenum", "Mo", "42")

            atom = atom_t(name="molybdenum", symbol="Mo",    atomic_mass=95.904678,  znuclear=42, core=36, nvalence=6)

        case("technetium", "Technetium", "Tc", "43")

            atom = atom_t(name="technetium", symbol="Tc",    atomic_mass=96.906364,  znuclear=43, core=36, nvalence=7)

        case("ruthenium", "Ruthenium", "Ru", "44")

            atom = atom_t(name="ruthenium", symbol="Ru",    atomic_mass=101.904348, znuclear=44, core=36, nvalence=8)

        case("rhodium", "Rhodium", "Rh", "45")

            atom = atom_t(name="rhodium",   symbol="Rh",    atomic_mass=102.905500, znuclear=45, core=36, nvalence=9)

        case("palladium", "Palladium", "Pd", "46")

            atom = atom_t(name="palladium", symbol="Pd",    atomic_mass=105.903478, znuclear=46, core=36, nvalence=10)

        case("silver", "Silver", "Ag", "47")

            atom = atom_t(name="silver",    symbol="Ag",    atomic_mass=106.905092, znuclear=47, core=36, nvalence=11)

        case("cadmium", "Cadmium", "Cd", "48")

            atom = atom_t(name="cadmium",   symbol="Cd",    atomic_mass=113.903357, znuclear=48, core=36, nvalence=12)

        case("indium", "Indium", "In", "49")

            atom = atom_t(name="indium",    symbol="In",    atomic_mass=114.903800, znuclear=49, core=36, nvalence=13)

        case("tin", "Tin", "Sn", "50")

            atom = atom_t(name="tin",       symbol="Sn",    atomic_mass=119.902220, znuclear=50, core=36, nvalence=14)

        case("antimony", "Antimony", "Sb", "51")

            atom = atom_t(name="antimony",  symbol="Sb",    atomic_mass=120.903821, znuclear=51, core=36, nvalence=15)

        case("tellurium", "Tellurium", "Te", "52")

            atom = atom_t(name="tellurium", symbol="Te",    atomic_mass=129.906229, znuclear=52, core=36, nvalence=16)

        case("iodine", "Iodine", "I", "53")

            atom = atom_t(name="iodine",    symbol="I",     atomic_mass=126.904473, znuclear=53, core=36, nvalence=17)

        case("xenon", "Xenon", "Xe", "54")

            atom = atom_t(name="xenon",     symbol="Xe",    atomic_mass=131.904144, znuclear=54, core=36, nvalence=18)

        case("cesium", "Cesium", "Cs", "55")

            atom = atom_t(name="cesium",    symbol="Cs",    atomic_mass=132.905429, znuclear=55, core=54, nvalence=1)

        case("barium", "Barium", "Ba", "56")

            atom = atom_t(name="barium",    symbol="Ba",    atomic_mass=137.905232, znuclear=56, core=54, nvalence=2)

        case("lanthanum", "Lanthanum", "La", "57")

            atom = atom_t(name="lanthanum", symbol="La",    atomic_mass=138.906346, znuclear=57, core=54, nvalence=3)

        case("cerium", "Cerium", "Ce", "58")

            atom = atom_t(name="cerium",    symbol="Ce",    atomic_mass=139.905433, znuclear=58, core=54, nvalence=4)

        case("praseodymium", "Praseodymium", "Pr", "59")

            atom = atom_t(name="praseodymium", symbol="Pr", atomic_mass=140.907647, znuclear=59, core=54, nvalence=5)

        case("neodymium", "Neodymium", "Nd", "60")

            atom = atom_t(name="neodymium", symbol="Nd",    atomic_mass=141.907719, znuclear=60, core=54, nvalence=6)

        case("promethium", "Promethium", "Pm", "61")

            atom = atom_t(name="promethium", symbol="Pm",   atomic_mass=144.912744, znuclear=61, core=54, nvalence=7)

        case("samarium", "Samarium", "Sm", "62")

            atom = atom_t(name="samarium",  symbol="Sm",    atomic_mass=151.919729, znuclear=62, core=54, nvalence=8)

        case("europium", "Europium", "Eu", "63")

            atom = atom_t(name="europium",  symbol="Eu",    atomic_mass=152.921225, znuclear=63, core=54, nvalence=9)

        case("gadolinium", "Gadolinium", "Gd", "64")

            atom = atom_t(name="gadolinium", symbol="Gd",   atomic_mass=157.924099, znuclear=64, core=54, nvalence=10)

        case("terbium", "Terbium", "Tb", "65")

            atom = atom_t(name="terbium",   symbol="Tb",    atomic_mass=158.925342, znuclear=65, core=54, nvalence=11)

        case("dysprosium", "Dysprosium", "Dy", "66")

            atom = atom_t(name="dysprosium", symbol="Dy",   atomic_mass=163.929171, znuclear=66, core=54, nvalence=12)

        case("holmium", "Holmium", "Ho", "67")

            atom = atom_t(name="holmium",   symbol="Ho",    atomic_mass=164.93032, znuclear=67, core=54, nvalence=13)

        case("erbium", "Erbium", "Er", "68")

            atom = atom_t(name="erbium",    symbol="Er",    atomic_mass=165.930290, znuclear=68, core=54, nvalence=14)

        case("thulium", "Thulium", "Tm", "69")

            atom = atom_t(name="thulium",   symbol="Tm",    atomic_mass=168.93421, znuclear=69, core=54, nvalence=15)

        case("ytterbium", "Ytterbium", "Yb", "70")

            atom = atom_t(name="ytterbium", symbol="Yb",    atomic_mass=173.938859, znuclear=70, core=54, nvalence=16)

        case("lutetium", "Lutetium", "Lu", "71")

            atom = atom_t(name="lutetium",  symbol="Lu",    atomic_mass=174.940770, znuclear=71, core=54, nvalence=17)

        case("hafnium", "Hafnium", "Hf", "72")

            atom = atom_t(name="hafnium",   symbol="Hf",    atomic_mass=179.946545, znuclear=72, core=54, nvalence=18)

        case("tantalum", "Tantalum", "Ta", "73")

            atom = atom_t(name="tantalum",  symbol="Ta",    atomic_mass=220.011368, znuclear=73, core=54, nvalence=19)

        case("tungsten", "Tungsten", "W", "74")

            atom = atom_t(name="tungsten",  symbol="W",     atomic_mass=183.95093, znuclear=74, core=54, nvalence=20)

        case("rhenium", "Rhenium", "Re", "75")

            atom = atom_t(name="rhenium",   symbol="Re",    atomic_mass=186.955765, znuclear=75, core=54, nvalence=21)

        case("osmium", "Osmium", "Os", "76")

            atom = atom_t(name="osmium",    symbol="Os",    atomic_mass=191.961479, znuclear=76, core=54, nvalence=22)

        case("iridium", "Iridium", "Ir", "77")

            atom = atom_t(name="iridium",   symbol="Ir",    atomic_mass=192.962924, znuclear=77, core=54, nvalence=23)

        case("platinum", "Platinum", "Pt", "78")

            atom = atom_t(name="platinum",  symbol="Pt",    atomic_mass=195.965833, znuclear=78, core=54, nvalence=24)

        ! The number of core and valence electrons for gold has beed changed for the tests. These values to be overwritten by
        ! the vallues from the ECP or pseudopotential files.
        case("gold", "Gold", "Au", "79")

            atom = atom_t(name="gold",      symbol="Au",    atomic_mass=196.966543, znuclear=79, core=60, nvalence=19)

        case("mercury", "Mercury", "Hg", "80")

            atom = atom_t(name="mercury",   symbol="Hg",    atomic_mass=200.967442, znuclear=80, core=54, nvalence=26)

        case("thallium", "Thallium", "Tl", "81")

            atom = atom_t(name="thallium",  symbol="Tl",    atomic_mass=204.974412, znuclear=81, core=54, nvalence=27)

        case("lead", "Lead", "Pb", "82")

            atom = atom_t(name="lead",      symbol="Pb",    atomic_mass=207.976636, znuclear=82, core=54, nvalence=28)

        case("bismuth", "Bismuth", "Bi", "83")

            atom = atom_t(name="bismuth",   symbol="Bi",    atomic_mass=208.980383, znuclear=83, core=54, nvalence=29)

        case("polonium", "Polonium", "Po", "84")

            atom = atom_t(name="polonium",  symbol="Po",    atomic_mass=209.987131, znuclear=84, core=54, nvalence=30)

        case("astatine", "Astatine", "At", "85")

            atom = atom_t(name="astatine",  symbol="At",    atomic_mass=210.990655, znuclear=85, core=54, nvalence=31)

        case("radon", "Radon", "Rn", "86")

            atom = atom_t(name="radon",     symbol="Rn",    atomic_mass=222.017570, znuclear=86, core=54, nvalence=32)

        case default
                error stop "Unknown element's atomic number or symbol. Atomic number > 86 (elements after Radon) not yet supported."
    end select
    return
    end function

end module
