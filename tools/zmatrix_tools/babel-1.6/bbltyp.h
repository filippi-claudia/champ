/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : bbl_types.h
AUTHOR(S) : Pat Walters
DATE : Constantly evolving

******/

#ifndef __BABEL_BBLTYP_H__
#define __BABEL_BBLTYP_H__

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "bblmacs.h"

enum file_type {
  none, mopac_cartesian, mopac_internal, mopac_output, macromodel, 
  mm2_output, pdb, csd_gstat, xyz, chemdraw, 
  idatm, mm2_input, alchemy, csd_cssr, free_fract,
  mac_molecule, report, amber_prep, molin, boogie, 
  caccrt, ball_and_stick, chem3d1, chem3d2, micro_world, 
  molfile, gamess_output, gamess_input, mmads, diagnostics, 
  shelx, csd_fdat, charmm, gaussian_output, 
  gaussian_input, gaussian_cart, hyperchem_hin, biosym_car, cadpac,
  cache_out, mm3, sybyl_mol, sybyl_mol2, spartan, quanta, molgen, wizard,
  feature, pcmodel, spart_semi, cacao_int, fenske_zmat,spart_mm,
  gaussian_template,bmin_com, conjure_tmplt, icon8, maccs, maccs2d,
  maccs3d, xed, unixyz, isis, dock, dock_pdb, bgf, m3d,
  psgvb_input, psgvb_output, psgvb_cart, psgvb_zmat,
  csr, smiles, gaussian_94, torlist, schakal, mol_inventor,
  gr96A,gr96N, tinker,box
  };

enum op_type {noaction, input, output};

enum type_err {zero, dummy, all_caps};

typedef struct
{
  enum file_type type;
  char *type_name;
  char *code;
  enum op_type operation;
#ifdef __cplusplus
  int (*func)(FILE *fp,struct ums *mol);
#else
  int (*func)();
#endif
} babel_rec;

typedef char char3[3];

typedef struct
{
  int a,b,c,d;
} torsion_rec;

typedef struct
{
  int a,b,c;
} angle_rec;

typedef struct
{
  int setlen;
  int *set;
} set_type;

typedef struct 
{
  int number;
  char name[3];
  double cov_rad;
  double bond_ord_rad;
  double vdw_rad;
  double bs_rad;
  int max_bonds;
  int color;
  double red;
  double green;
  double blue;
} element_type;

typedef struct 
{
  double x,y,z;
} coord_type;


typedef struct 
{
  coord_type point;
  char type[5];
  int max_bonds;
  int valence;
  int atomic_number;
  int connected_atoms[MAX_CONNECTIONS];
  int bond_order[MAX_CONNECTIONS];
  double radius;
  double bond_ord_rad;
  double dble;
  int organic;
  int redo;
  int pos[3];
  double charge;
  set_type *atm_set;
} atom_type;


typedef struct 
{
  int start;
  int end;
  int bond_order;
} connect_type;

typedef struct 
{
  double r;
  double w;
  double t;
/*  int  n; */
  int  na;
  int  nb;
  int  nc;
} int_type;

typedef struct
{
  char default_extension[BUFF_SIZE];
  char base_name[BUFF_SIZE];
  char infile_name[BUFF_SIZE];
  char outfile_name[BUFF_SIZE];
  char input_keywords[BUFF_SIZE];
  char output_keywords[BUFF_SIZE];
  char del_str[BUFF_SIZE];
  babel_rec input_info;
  babel_rec output_info;
  int verbose;
  int use_menus;
  int do_add_hydrogens;
  int do_delete_atoms;
  int the_size;
  int multi;
  int lower_limit;
  int upper_limit;
  int spline;
  int increment;
  int renum;
  int new_base;
  int center;
  int align;
  int precipitate;
  int new_file;
  int calc_charges;
  int orient;
  int no_dummy;
  int push_hydrogens;
} bbl_control;



/*---------- structures for ring detection --- */

typedef struct
{
  int level;
  int ancestor;
} spanning_tree;

typedef struct
{
  int length;
  int bogus;
  set_type *path_set;
  int *path_atoms;
  int closure;
  int found;
} path;

typedef struct 
{
  int count;
  path *ring_list;
} ring_struct;


typedef struct
{
  int num;
  set_type *ring_atms;
  set_type *arom_atms;
  set_type *arom_rings;
  set_type **rings;
} ring_info;

/*--------------------------------------------*/


/*************
HUGEPTR as defined in the struct below is a maco
which defines HUGEPTR as the keyword huge when #define MSDOS
is present or as nothing when #define MSDOS is absent.  This is 
necessary because of the lovely method Intel uses to segment memory.

The macros are defined in bblmacs.h
PW - 072094
****************/

typedef struct
{
  double A,B,C;
  double Alpha,Beta,Gamma;
} fract_type;

typedef struct
{
  int serial_num;
  int chain_num;
  int res_num;
  char res_type[6];
  char atm_type[6];
} res_list;

struct ums
{
  bbl_control *control;
  ring_struct *rings;
  ring_info *ring_info;
  char default_extension[20];
  int num_atoms;
  int num_bonds;
  int num_connections_initialized;
  int num_residues_initialized;
  int num_internal_initialized;
  char title[BUFF_SIZE];
  double energy;
  res_list HUGEPTR *residues;
  atom_type HUGEPTR *atoms;
  connect_type HUGEPTR *connections;
  int_type HUGEPTR *internal;
  fract_type HUGEPTR *fract;
  struct ums *next;
};
typedef struct ums ums_type;

/* structs for miniums functions */

struct mini_ums_type
{
  int num_atoms;
  double energy;
  coord_type *atoms;
  struct mini_ums_type *next;
};

typedef struct mini_ums_type mini_ums;

typedef struct 
{
  int atom_number;
  int heavy_atoms;
} heavy_rec;

typedef struct 
{
  char name[8];
  int number;
} pdb_type_rec;


typedef struct
{
  char *text;
  enum file_type type;
  int keywords_required;
} menu_type;


typedef struct
{
  double x,y,z;
} vect_type;


typedef struct
{
  double a1, b1, c1;
  double a2, b2, c2;
  double a3, b3, c3;
} matrix_3x3;


typedef struct
{
  double x;
  double y;
  double z;
  int num;
  double dist;
} temp_atom_rec;

typedef char warning[80];

typedef struct 
{
  char sym[10];
  double val;
} zsymbol;

typedef struct
{
  int level;
  int ancestor;
  int kids;
  int kid[MAX_CONNECTIONS];
} z_tree;


/* structure and macros for set operations */
#define SETWORD  32  /* platform dependant - should == sizeof(int) */

#define LowBit(set, bit)\
  {register int m;\
   if (set != 0)\
   {\
      bit = 31;\
      if (set != 0x80000000) {\
      if (m = (set & 0x0000ffff)) {set = m; bit -= 16;}\
      if (m = (set & 0x00ff00ff)) {set = m; bit -= 8;}\
      if (m = (set & 0x0f0f0f0f)) {set = m; bit -= 4;}\
      if (m = (set & 0x33333333)) {set = m; bit -= 2;}\
      if (m = (set & 0x55555555)) {set = m; bit -= 1;}}}\
   else bit = -1;}

static int bitsoff[SETWORD] =
{
0xFFFFFFFF,0xFFFFFFFE,0xFFFFFFFC,0xFFFFFFF8,0xFFFFFFF0,0xFFFFFFE0,0xFFFFFFC0,
0xFFFFFF80,0xFFFFFF00,0xFFFFFE00,0xFFFFFC00,0xFFFFF800,0xFFFFF000,0xFFFFE000,
0xFFFFC000,0xFFFF8000,0xFFFF0000,0xFFFE0000,0xFFFC0000,0xFFF80000,0xFFF00000,
0xFFE00000,0xFFC00000,0xFF800000,0xFF000000,0xFE000000,0xFC000000,0xF8000000,
0xF0000000,0xE0000000,0xC0000000,0x80000000
};


#define biton(seta,member)   seta->set[(member / SETWORD)] |= \
(1 << (member % SETWORD));

#define bitoff(seta,member)   seta->set[(member / SETWORD)] &= \
(~(1 << (member % SETWORD)));

#define bit_is_on(seta,member) (seta->set[(member/SETWORD)]>>\
                                (member % SETWORD)&1)

enum coord_state {linear, trigonal, tetrahedral, octahedral};

typedef struct
{
  int ptr;
  int *bond;
  int *choice;
} bnd_stack;

typedef struct
{
  double is,id,it,wsd,wdt;
} bnd_info;

typedef struct 
{
  char atm_typ[3];
  double a,b,c,d;
} gast_param;



/*----------------------------------------------------
Data structures for Simon Kilvington's SMILES routines
-----------------------------------------------------*/

typedef void *block_ptr;

/* smiles connection table data */

#define md_MAXBONDS  8	    /* max atoms that an atom can be bonded to */
#define md_NOBOND   -1	    /* used in bondedto[] arrays */

typedef char smilesbond_t;

#define SMILESBOND_NoBond   ((smilesbond_t) 0)
#define SMILESBOND_Single   ((smilesbond_t) 1)
#define SMILESBOND_Double   ((smilesbond_t) 2)
#define SMILESBOND_Triple   ((smilesbond_t) 3)
#define SMILESBOND_Aromatic ((smilesbond_t) 5)

typedef struct
{
   char         symbol[2];
   int          bondedto[md_MAXBONDS];	    /* indices of atoms this one is bonded to */
   smilesbond_t bondtype[md_MAXBONDS];
} smilesatom_t;

typedef struct
{
   int          natoms;
   smilesatom_t atom[1];		    /* array extends to atom[natoms-1] */
} smilescontab_t;

#define smiles_CONTABHDRSIZE (sizeof(smilescontab_t) - sizeof(smilesatom_t))

/*-------------------------------------------------------------------------------------*/


#include "babel.h"

#endif /* !__BABEL_BBLTYP_H__ */
