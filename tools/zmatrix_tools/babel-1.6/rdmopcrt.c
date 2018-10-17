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

FILE : rdmopcrt.c
AUTHOR(S) : Pat Walters
DATE : 10-92    modified (very slightly!) 7-94
PURPOSE : Routines to read a mopac cartesian coordinate input file

******/

#include "bbltyp.h"

int 
read_mop_cart(FILE *file1, ums_type *mol)
{
  char mop_cart_line[BUFF_SIZE];
  int count = 0;
  int i;
  int result;
  
  for (i = 0; i < 3; i ++)
  {
    fgets(mop_cart_line,sizeof(mop_cart_line), file1);
    if ( i == 2)
      strcpy(Title, mop_cart_line);
  }
  while (fgets(mop_cart_line,sizeof(mop_cart_line), file1) != NULL)
  {
    if (count_tokens(mop_cart_line," \t\n") >= 7)
	{
	  count++;
	}
  }
  Atoms = count;
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  fseek(file1,0,0);
  count = MIN_ATOM;
  for (i = 0; i < 3; i ++)
    fgets(mop_cart_line,sizeof(mop_cart_line), file1);
  while (fgets(mop_cart_line,sizeof(mop_cart_line), file1) != NULL)
  {
    if (count_tokens(mop_cart_line," \t\n") >= 7) 
    {
      UpdateProgress();
      sscanf(mop_cart_line,"%s%lf%d%lf%d%lf%d",
	     Type(count),
	     &X(count),
	     &mol->atoms[count].pos[0],
	     &Y(count),
	     &mol->atoms[count].pos[1],
	     &Z(count),
	     &mol->atoms[count].pos[2]);
      clean_atom_type(Type(count));      
      count ++;
    }
  }
  mol->atoms[0].pos[0] = 0;

  if (Atoms > 0)
  {
    result = assign_radii(mol);
    result = assign_bonds(mol);
    result = assign_types(mol);
    assign_bond_order(mol);
  }
  return(TRUE);  
}

   
    
    
    
	  


