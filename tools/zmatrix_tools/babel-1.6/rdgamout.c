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

FILE : rgamout.c
AUTHOR(S) : Pat Walters
DATE : 12-93
PURPOSE : Routines to read a GAMESS output file

******/

#include "bbltyp.h"

extern element_type *elements;
#define BOHR_TO_ANGSTROM 0.5292

int 
read_gamess_output(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  char temp_type[5];
  int num;
  char xstr[15], ystr[15], zstr[15];
  int i = 0;
  int result;
  long pos = 0;
  int tokens;
  int optimized = FALSE;

  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (strstr(the_line,"COORDINATES OF ALL ATOMS ARE") != NULL)
    {
      optimized = TRUE;
      show_warning("Using geometry optimized coordinates");
      pos = ftell(file1);
      fgets(the_line,sizeof(the_line), file1);
      fgets(the_line,sizeof(the_line), file1);
      Atoms = -1;
      tokens = 5;
      while (tokens == 5) 
      {
	if (fgets(the_line,sizeof(the_line), file1) != NULL)
	  tokens = -1;
	tokens = count_tokens(the_line," \t\n");
	Atoms++;
      }
    }
  }
  
  if (!optimized)
  {
    tokens = 5;
    rewind(file1);
    while (fgets(the_line,sizeof(the_line), file1) != NULL)
    {
      if (strstr(the_line,"ATOMIC                      COORDINATES (BOHR)") != NULL)
      {
	Atoms = 0;
	pos = ftell(file1);
	show_warning("Could not find optimized coordinates");
	show_warning("Using input coordinates");
	fgets(the_line,sizeof(the_line), file1);
	fgets(the_line,sizeof(the_line), file1);
	while (tokens == 5) 
	{
	  if (fgets(the_line,sizeof(the_line), file1) != NULL)
	    tokens = -1;
	  tokens = count_tokens(the_line," \t\n");
	  Atoms++;
	}
      }
    }
  }
  
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  fseek(file1,pos,0);
  fgets(the_line,sizeof(the_line), file1);
  if (optimized)
    fgets(the_line,sizeof(the_line), file1);
  for (i = 1; i <= Atoms; i++)
  {
    UpdateProgress();
    fgets(the_line,sizeof(the_line), file1);
    sscanf(the_line,"%*s%s%s%s%s",temp_type,xstr,ystr,zstr); 
    X(i) = (double)atof(xstr);
    Y(i) = (double)atof(ystr);
    Z(i) = (double)atof(zstr);
    num = atoi(temp_type);
    strcpy(Type(i),elements[num].name);
  }
  if (!optimized)
    bohr_to_angstroms(mol);
  
  if (Atoms > 0)
  {
    result = assign_radii(mol);
    result = assign_bonds(mol);
    result = assign_types(mol);
    assign_bond_order(mol);
    result = build_connection_table(mol);
    assign_bond_order(mol);
  }
  read_to_eof(file1);

  return(TRUE);
}

       
void bohr_to_angstroms(ums_type *mol)
{
  int i;
  
  for (i = 1; i <= Atoms; i++)
  {
    X(i) = X(i) * BOHR_TO_ANGSTROM;
    Y(i) = Y(i) * BOHR_TO_ANGSTROM;
    Z(i) = Z(i) * BOHR_TO_ANGSTROM;
  }
}























