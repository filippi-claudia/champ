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

FILE : rdgauout.c
AUTHOR(S) : Pat Walters
DATE : 2-20-96
PURPOSE : routines to read a Gaussian 94 log file
******/

#include "bbltyp.h"

#define SEARCH_STRING "Standard orientation:"

int read_gaussian_94(FILE *fp, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  long pos = 0;
  int i;
  
  Atoms = 0;
  while (fgets(the_line,sizeof(the_line), fp) != NULL)
  {
    if (strstr(the_line,SEARCH_STRING))
      break;
  }
  toss(fp,4);  /* throw away four lines */
  pos = ftell(fp);  /* see where we are */

  /* find out how many atoms there are */
  while (fgets(the_line,sizeof(the_line), fp) != NULL)
  {
    if (strstr(the_line,"-----------------"))
      break;
    Atoms++;
  }
  
  initialize_ums(&mol);
  fseek(fp,pos,0);
  for (i = 1; i <= Atoms; i++)
  {
    fgets(the_line,sizeof(the_line), fp);
    sscanf(the_line,"%*s%d%lf%lf%lf",&Atomic_number(i),&X(i),&Y(i),&Z(i));
  }
  if (Atoms > 0)
  {
    add_element_types(mol);
    assign_radii(mol);
    assign_bonds(mol);
    assign_types(mol);
    assign_bond_order(mol);
  }

  pos = ftell(fp);

  while (fgets(the_line,sizeof(the_line), fp) != NULL) 
  {
    if (strstr(the_line,SEARCH_STRING))
    {
      fseek(fp,pos,0);
      break;
    }
  }
  
  return(TRUE);
}

#undef SEARCH_STRING  
                   
