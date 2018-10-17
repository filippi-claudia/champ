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

FILE : rdprep.c
AUTHOR(S) : Pat Walters
DATE : 11-93
PURPOSE : Routines to read an AMBER PREP file

******/

#include "bbltyp.h"

int 
read_amber_prep(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i = 0;
  int result;
  ums_type *new_mol;

  new_mol = (ums_type *)malloc(sizeof(ums_type));
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (count_tokens(the_line," ") > 8)
    {
      i++;
    }
  }
  Atoms = i;
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  initialize_internal(&mol);
  rewind(file1);
  i = 1;
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (count_tokens(the_line," ") > 8)
    {
      UpdateProgress();
      sscanf(the_line,"%*d %s %*s %*s %d %d %d %lf %lf %lf %*f",
	     Type(i),
	     &mol->internal[i].na,
	     &mol->internal[i].nb,
	     &mol->internal[i].nc,
	     &mol->internal[i].r,
	     &mol->internal[i].w,
	     &mol->internal[i].t);
      if (strchr("C H O N S",Type(i)[0]) != NULL)
      {
	Type(i)[1] = '\0';
      }
      clean_atom_type(Type(i));
      if EQ(Type(i),"Du")
	strcpy(Type(i),"DU");
      i++;
    }
  }
  int_to_cart(mol);
  result = assign_radii(mol);
  result = assign_bonds(mol);
  result = assign_types(mol);
  result = build_connection_table(mol);
  assign_bond_order(mol);

  new_mol = delete_atoms(mol,"Du");
  Atoms = new_mol->num_atoms;
  Bonds = new_mol->num_bonds;
  mol->atoms = new_mol->atoms;
  mol->connections = new_mol->connections;
  return(TRUE);  
}

   





    
    
	  










