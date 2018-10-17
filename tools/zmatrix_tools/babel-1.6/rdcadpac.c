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

FILE : rdcadpac.c
AUTHOR(S) : Pat Walters
DATE : 2-94
PURPOSE : Routines to read a CADPAC output file

******/

#include "bbltyp.h"

int 
  read_cadpac(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  long pos = 0;
  int optimized = FALSE;
  int tokens = 5;
  int i;
  int result;
  
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (strstr(the_line,"New Molecular Geometry") != NULL)
    {
      optimized = TRUE;
      pos = ftell(file1);
    }
  }
#ifdef MAC
ShowProgress(Atoms,"Reading Atoms");
#endif
  if (optimized)
  {
    rewind(file1);
    fseek(file1,pos,0);
    fgets(the_line,sizeof(the_line),file1);
    fgets(the_line,sizeof(the_line),file1);
    
    Atoms = -1;
    while (tokens == 5)
    {
      fgets(the_line,sizeof(the_line),file1);
      tokens = count_tokens(the_line,"\t\n ");
      Atoms++;
    }
    initialize_ums(&mol);
    rewind(file1);
    fseek(file1,pos,0);
    fgets(the_line,sizeof(the_line),file1);
    fgets(the_line,sizeof(the_line),file1);
    for (i = 1; i <= Atoms; i++)
    {
#ifdef MAC
      UpdateProgress();
#endif
      fgets(the_line,sizeof(the_line),file1);
      sscanf(the_line,"%*d%s%lf%lf%lf",
	     Type(i),
	     &X(i),
	     &Y(i),
	     &Z(i));
      clean_atom_type(Type(i));
    }  
  }

  else
  {
    rewind(file1);
    strcpy(the_line,"");
    while (strstr(the_line,"Geometry, (in atomic units)") == NULL)
    {
      fgets(the_line,sizeof(the_line),file1);
    }
    for (i = 0; i < 5; i++)
      fgets(the_line,sizeof(the_line),file1);
    pos = ftell(file1);
    Atoms = -1;
    while (tokens == 5)
    {
      fgets(the_line,sizeof(the_line),file1);
      tokens = count_tokens(the_line,"\t\n ");
      Atoms++;
    }
    initialize_ums(&mol);
    fseek(file1,pos,0);
    for (i = 1; i <= Atoms; i++)
    {
#ifdef MAC
      UpdateProgress();
#endif
      fgets(the_line,sizeof(the_line),file1);
      sscanf(the_line,"%s%*f%lf%lf%lf",
	     Type(i),
	     &X(i),
	     &Y(i),
	     &Z(i));
      clean_atom_type(Type(i));
    }  
  }

  bohr_to_angstroms(mol);
  if (Atoms > 0)
  {
    result = assign_radii(mol);
    result = assign_bonds(mol);
    result = assign_types(mol);
    result = build_connection_table(mol);
    assign_bond_order(mol);
  }
  return(TRUE);
}





