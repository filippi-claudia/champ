/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------
FILE : bbl_macros.h
AUTHOR(S) : Pat Walters
DATE : 8-93
PURPOSE : Contains macros used by the babel program

******/

#ifndef __BABEL_BBLMACS_H__
#define __BABEL_BBLMACS_H__

#define BABEL_VERSION "1.6 "
#define STARS "************************"
#define NOKEY "KEYWORDS GO HERE"

#ifdef MSDOS
#define HUGEPTR huge
#else
#define HUGEPTR 
#endif


#define NULL_CHAR '\0'
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define BUFF_SIZE 300
#define NUM_TYPES 31
#define MAX_CONNECTIONS 20
#define MAX_ELEMENTS 108
#define MIN_ATOM 1

#ifndef PI
#define PI 3.1415926535897932384626433
#endif

#define RAD_TO_DEG 180.0/PI
#define DEG_TO_RAD PI/180.0

#ifndef ULTRIX
#define QSORT_PROTO (int(*)(const void *,const void *))
#else
#define QSORT_PROTO
#endif

#define MASTERSIZE sizeof(master) / sizeof(babel_rec)

/* Color definitions */

#define BBL_UNDEF   0
#define BBL_BLACK   1
#define BBL_GREY    2
#define BBL_DKBLU   3
#define BBL_BLUE    4
#define BBL_LTBLU   5
#define BBL_AQUA    6
#define BBL_TURQ    7
#define BBL_BLUGRN  8
#define BBL_DKGRN   9
#define BBL_GREEN   10
#define BBL_LTGRN   11
#define BBL_YELGRN  12
#define BBL_YELLOW  13
#define BBL_ORANGE  14
#define BBL_DKRED   15
#define BBL_RED     16
#define BBL_PINK    17
#define BBL_REDPUR  18
#define BBL_PURPLE  19
#define BBL_BLUPUR  20
#define BBL_WHITE   21

/* Macro definitions */
#define SQUARE(x) (x) * (x)
#define NEW(type) (type *) malloc(sizeof(type))
#define EQ(a, b)        (strcmp((a), (b)) == 0)
#define NOTEQ(a, b)     (strcmp((a), (b)) != 0)
#define EQn(a, b, n)    (strncmp((a), (b), (n)) == 0)
#define NOTEQn(a, b, n)    (strncmp((a), (b), (n)) != 0)

#ifndef MIN
#define MIN(a,b) (((a) < (b) ? (a) : (b)))
#endif
#ifndef MAX
#define MAX(a,b) (((a) > (b) ? (a) : (b)))
#endif

#define SWAP(a,b) {a ^= b; b ^= a; a ^= b;}

/* Macros defining UMS components */
#define Type(x)             mol->atoms[x].type
#define Valence(x)          mol->atoms[x].valence
#define Max_bonds(x)        mol->atoms[x].max_bonds
#define Connection(x,y)     mol->atoms[x].connected_atoms[y]
#define Atoms               mol->num_atoms
#define Point(x)            mol->atoms[x].point
#define Start(x)            mol->connections[x].start
#define End(x)              mol->connections[x].end
#define Bond_order(x)       mol->connections[x].bond_order
#define Bonds               mol->num_bonds
#define X(a)                mol->atoms[a].point.x
#define Y(a)                mol->atoms[a].point.y
#define Z(a)                mol->atoms[a].point.z
#define BO(x,y)             mol->atoms[x].bond_order[y]
#define Redo(x)             mol->atoms[x].redo
#define Radius(x)           mol->atoms[x].radius
#define Atomic_number(x)    mol->atoms[x].atomic_number
#define BORadius(x)         mol->atoms[x].bond_ord_rad
#define Energy              mol->energy
#define Double(x)           mol->atoms[x].dble
#define Next                mol->next
#define Title               mol->title
#define Charge(x)           mol->atoms[x].charge
#define Organic(x)          mol->atoms[x].organic
#define IsOrganic(x)        (mol->atoms[x].organic == 1)

/*Macros defining the control structure */
#define DefaultExtension    mol->control->default_extension

#define InfileName          mol->control->infile_name
#define InputKeywords       mol->control->input_keywords
#define InputInfo           mol->control->input_info
#define InfileType          mol->control->input_info.type
#define InputTypeName       mol->control->input_info.type_name
#define ReaderFunction      mol->control->input_info.func
#define InputTrans          mol->control->input_info.translate

#define OutfileName         mol->control->outfile_name
#define OutputKeywords      mol->control->output_keywords
#define OutputInfo          mol->control->output_info
#define OutfileType         mol->control->output_info.type
#define OutfileTypeName     mol->control->output_info.type_name
#define WriterFunction      mol->control->output_info.func
#define OutputTrans         mol->control->output_info.translate

#define BaseName            mol->control->base_name
#define Verbose             mol->control->verbose
#define AddHydrogens        mol->control->do_add_hydrogens
#define DeleteAtoms         mol->control->do_delete_atoms
#define UseMenus            mol->control->use_menus
#define Size                mol->control->the_size
#define Multi               mol->control->multi
#define DeleteStr           mol->control->del_str
#define LowerLimit          mol->control->lower_limit
#define UpperLimit          mol->control->upper_limit
#define Spline              mol->control->spline
#define Increment           mol->control->increment
#define DoRenum             mol->control->renum
#define NewBase             mol->control->new_base
#define CenterMol           mol->control->center
#define Align               mol->control->align
#define StdOrientation      mol->control->orient
#define Precipitate         mol->control->precipitate
#define PushHydrogens       mol->control->push_hydrogens
#define MakeNewFile         mol->control->new_file
#define CalcCharges         mol->control->calc_charges
#define NoDummy             mol->control->no_dummy

#define NA(x) mol->internal[x].na
#define NB(x) mol->internal[x].nb
#define NC(x) mol->internal[x].nc
#define R(x) mol->internal[x].r
#define W(x) mol->internal[x].w
#define T(x) mol->internal[x].t

#define SerialNum(x) mol->residues[x].serial_num
#define ChainNum(x) mol->residues[x].chain_num
#define ResNum(x) mol->residues[x].res_num
#define ResName(x) mol->residues[x].res_type
#define AtmId(x) mol->residues[x].atm_type
#define HasResidues (mol->residues != NULL)

#define single_struct         0
#define multi_struct          1
#define multi_conf            2
#define sequential_name       3
#define title_as_name         4

#define SP3_MAX      114.0  /*changed from 115.0 PW->5-17-93*/
#define MAY_BE_SP2   122.0
#define SP_MIN       160.0

#define V1_C1_C1_CUTOFF 1.22
#define V1_C2_C_CUTOFF  1.41
#define V1_C2_N_CUTOFF  1.37

#define V1_N1_C1_CUTOFF 1.20 
#define V1_N3_C_CUTOFF  1.38
#define V1_N3_N3_CUTOFF 1.43
#define V1_N3_N2_CUTOFF 1.41

#define V1_O2_C2_CUTOFF 1.30
#define V1_O2_AS_CUTOFF 1.685

#define V1_S2_C2_CUTOFF 1.76
#define V1_S2_AS_CUTOFF 2.11

#define V2_C3_C_CUTOFF  1.53
#define V2_C3_N_CUTOFF  1.46
#define V2_C3_O_CUTOFF  1.44

#define V2_N2_C_CUTOFF  1.38
#define V2_N2_N_CUTOFF  1.32

#define V2_C2_C_CUTOFF  1.42
#define V2_C2_N_CUTOFF  1.41
#define GEN_C3_C_CUTOFF 1.45

#define RingSize(x)   rings->ring_list[x].length
#define NumRings      rings->count
#define RingAtom(x,y) rings->ring_list[x].path_atoms[y]

#define mX(a) mini->atoms[a].x
#define mY(a) mini->atoms[a].y
#define mZ(a) mini->atoms[a].z
#define mAtoms mini->num_atoms
#define mPoint(a) mini->atoms[a]
#define mEnergy mini->energy

#define ONE_OVER_SQRT3  0.577350269
#define SQRT_TWO_THIRDS 0.816496581

#define SP3_C_H_DIST 1.115
#define SP2_C_H_DIST 1.103
#define SP_C_H_DIST 1.090

#define SP3_N_H_DIST 1.020
#define SP2_N_H_DIST 1.020

#define SP3_O_H_DIST 0.950

#define IsUnsatType(x)  (EQ(x,"Car") || EQ(x,"C2") || EQ(x,"Sox") || EQ(x,"Sac") || EQ(x,"Pac") || EQ(x,"So2"))

#endif  /* !__BABEL_BBLMACS_H__ */
