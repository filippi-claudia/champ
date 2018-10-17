/*****
  This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
  
  For more information please contact :
  
  babel@mercury.aichem.arizona.edu
  ---------------------------------------------------------------------------
  
  FILE : rdmopac.c
  AUTHOR(S) : Pat Walters
  DATE : 1-93
  PURPOSE : Routines to read a Spartan Semiempirical file
  ******/

#include "bbltyp.h"

int 
  read_spartan_semiempirical(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i = 0;
  int result;
  int ret = TRUE;
  long pos = 0;
  int tokens = 5;
  int done = FALSE;
  
  while ((fgets(the_line,sizeof(the_line), file1) != NULL) && (i < 2))
  {
    if (strstr(the_line,"Cartesian"))
      i++;
  }
  if (i > 1)
  {
    for (i = 0; i < 2; i++)
    {
      fgets(the_line,sizeof(the_line), file1);
    }
    pos = ftell(file1);
    Atoms = 0;
    while (tokens == 5)
    {
      fgets(the_line,sizeof(the_line), file1);
      tokens = count_tokens(the_line," \n\t");
      if (tokens == 5)
	Atoms++;
    }
    initialize_ums(&mol);
    fseek(file1,pos,0);
    for (i = 1; i <= Atoms; i++)
    {
      fgets(the_line,sizeof(the_line), file1);
      sscanf(the_line,"%s %*d %lf %lf %lf",Type(i),&X(i),&Y(i),&Z(i)); 
    }
    while ((fgets(the_line,sizeof(the_line), file1) != NULL) && (!done))
    {
      if (strstr(the_line,"Heat"))
      {
	sscanf(&the_line[23],"%lf",&Energy);
	done = TRUE;
      }
    }      
    if (Atoms > 0)
    {
      result = assign_radii(mol);
      result = assign_bonds(mol);
      result = assign_types(mol);
      result = build_connection_table(mol);
      assign_bond_order(mol);
    }
  }
  else
  {
    show_warning("No optimized coordinates in this file");
    ret = FALSE;
  }
  fclose(file1);
  return(ret);
}

       
























