! Enter licence information here


! Licence information ends here
!
! Note. The elemental data is taken from the NIST website (https://physics.nist.gov/PhysRefData/Handbook/periodictable_a.htm)
!       The data tabulated in this file is only for the most abundant isotope only.

module periodic_table
    public :: element, atom_t
    private :: i1, sp

    integer, parameter :: sp = kind(1.0)
    integer, parameter :: i1 = kind(int1)

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
                
        case default
                error stop "Unknown element's atomic number or symbol"
    end select
    return 
    end function      

end module