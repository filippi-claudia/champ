/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : rdunichem.c
AUTHOR(S) : Pat Walters
DATE : 2-28-95
PURPOSE : routines to read the Cray Unichem format

******/

#include "bbltyp.h"

int 
read_unichem(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i;
  int num;

  fgets(the_line,sizeof(the_line),file1);
  fgets(the_line,sizeof(the_line),file1);
  sscanf(the_line,"%d",&Atoms);
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);

  for (i = MIN_ATOM; i <= Atoms; i++)
  {
    UpdateProgress();
    fgets(the_line,sizeof(the_line), file1);
    sscanf(the_line,"%d %lf %lf %lf",
	   &num,
	   &X(i),
	   &Y(i),
	   &Z(i));
    atomic_number_to_name(num,Type(i));
    clean_atom_type(Type(i));
  }

  assign_radii(mol);
  assign_bonds(mol); 
  assign_types(mol);
  build_connection_table(mol); 
  assign_bond_order(mol);
    
  return(TRUE);
}








