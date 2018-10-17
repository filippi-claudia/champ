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

FILE : rdxyz.c
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : routines to read the XYZ format used by the Xmol program from MSC
	
MODIFIED : 10-16-93 to allow the use of multistructure files

******/


#include "bbltyp.h"

int 
read_xyz(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i;

  fgets(the_line,sizeof(the_line),file1);
  sscanf(the_line,"%d",&Atoms);
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  fgets(the_line,sizeof(the_line),file1);

  for (i = 1; i <= Atoms; i++)
  {
    UpdateProgress();
    fgets(the_line,sizeof(the_line), file1);
    sscanf(the_line,"%s %lf %lf %lf",
	   Type(i),
	   &X(i),
	   &Y(i),
	   &Z(i));
    clean_atom_type(Type(i));
  }

  assign_radii(mol);
  assign_bonds(mol); 
  assign_types(mol);
  build_connection_table(mol); 
  assign_bond_order(mol);
  return(TRUE);
}

