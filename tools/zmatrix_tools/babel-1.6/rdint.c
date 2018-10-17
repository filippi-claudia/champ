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

FILE : rdint.c
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : Routines to read a mopac internal coordinate file
******/

#include "bbltyp.h"

int 
read_mopint(FILE *file1, ums_type *mol)
{
  char mopint_line[BUFF_SIZE];
  int i = 0;
  int result;

/**** find out how many atoms are in the mopint file ****/

  for (i = 0; i < 3; i++)
    fgets(mopint_line,sizeof(mopint_line), file1);  
  i = 0;
  while (fgets(mopint_line,sizeof(mopint_line), file1) != NULL)
  {
    if (i < 3)
      i++;
    else
      if (count_tokens(mopint_line,"\t\n ") >= 10)
	i++;
  }
  Atoms = i;
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  initialize_internal(&mol);
  rewind(file1);
  for (i = 0; i < 3; i++)
    fgets(mopint_line,sizeof(mopint_line), file1);  
  i = MIN_ATOM;
  for (i = MIN_ATOM; i <= Atoms; i++)
  {
    fgets(mopint_line,sizeof(mopint_line), file1);  
    UpdateProgress();
    sscanf(mopint_line,"%s %lf %*d %lf %*d %lf %*d %d %d %d %*f",
	   Type(i),&R(i),&W(i),&T(i),&NA(i),&NB(i),&NC(i));
    clean_atom_type(Type(i));
  }

  R(1) = 0.0;
  W(1) = 0.0;
  T(1) = 0.0;
  NA(1) = 0;
  NB(1) = 0;
  NC(1) = 0;

  W(2) = 0.0;
  T(2) = 0.0;
  NB(2) = 0;
  NC(2) = 0;

  T(3) = 0.0;
  NC(3) = 0;

  if (Atoms > 0)
  {
    result = int_to_cart(mol);
    result = assign_radii(mol);
    result = assign_bonds(mol);
    result = assign_types(mol);
    result = build_connection_table(mol);
    assign_bond_order(mol);
  }
  read_to_eof(file1);
  return(TRUE);  
}

   
    
    
    
	  


