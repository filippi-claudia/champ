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
DATE : 5-96
PURPOSE : routine to read the file format used by the SCHAKAL program.
          This program was of course named for that famous basketball player
	  rapper and all around giant human Schakal O'Neill.
******/

#include "bbltyp.h"

int read_schakal(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i;
  long int pos;
  int tokens;
  matrix_3x3 m;
  fract_type f;
  int is_fractional = FALSE;
  
  fgets(the_line,sizeof(the_line),file1);
  sscanf(the_line,"%d",&Atoms);
  ShowProgress(Atoms,"Reading Atoms");
  
  pos = ftell(file1);
  while ((fgets(the_line,sizeof(the_line), file1) != NULL) && NOTEQn(the_line,"END",3))
  { 
    if (EQn(the_line,"CELL",4))
    {
      tokens = count_tokens(the_line,"\t \n");
      if (tokens == 7)
      {
	sscanf(the_line,"%*s%lf%lf%lf%lf%lf%lf",&f.A,&f.B,&f.C,&f.Alpha,&f.Beta,&f.Gamma);
	is_fractional = TRUE;
	fill_orth_matrix(&f,&m);
      }
    }
    if (EQn(the_line,"ATOM",4))
    {
      i++;
    }
  }
  
  fseek(file1,pos,0);
  Atoms = i;
  
  initialize_ums(&mol);
  i = 0;
  
  while ((fgets(the_line,sizeof(the_line), file1) != NULL) && NOTEQn(the_line,"END",3))
  { 
    if (EQn(the_line,"ATOM",4))
    {
      i++;
      sscanf(the_line,"%*s %s %lf %lf %lf",
	     Type(i),
	     &X(i),
	     &Y(i),
	     &Z(i));
      clean_atom_type(Type(i));
      if (is_fractional)
	fract_to_cart(&Point(i),&m);
    }
  }

  
  assign_radii(mol);
  assign_bonds(mol); 
  assign_types(mol);
  build_connection_table(mol); 
  assign_bond_order(mol);
  return(TRUE);
}

