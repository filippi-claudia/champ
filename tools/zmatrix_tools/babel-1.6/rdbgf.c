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
DATE : 12-2-95
PURPOSE : routines to read the BioDesign BGF format
******/

#include "bbltyp.h"

int 
read_bgf(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i,j;
  long pos;
  int con_atm, bond_order;
  int tokens;
  char temp_type[15];

  while (fgets(the_line,sizeof(the_line),file1))
  {
    if (EQn(the_line,"FORMAT",6))
      break;
  }

  ShowProgress(Atoms,"Reading Atoms");
  Atoms = 0;
  
  pos = ftell(file1);

  while (fgets(the_line,sizeof(the_line),file1))
  {
    if ((EQn(the_line,"ATOM",4)) || (EQn(the_line,"HETATM",6)))
      Atoms++;
    if (EQn(the_line,"FORMAT",6))
      break;
  }

  initialize_ums(&mol);
  fseek(file1,pos,0);

  for (i = 1; i <= Atoms; i++)
  {
    UpdateProgress();
    fgets(the_line,sizeof(the_line), file1);
    sscanf(the_line,"%*s %*s %*s %*s %*s %*s %lf %lf %lf %s %*s %*s %lf",
	   &X(i),
	   &Y(i),
	   &Z(i),
           temp_type,
	   &Charge(i));
    clean_atom_type(temp_type);
    strcpy(Type(i),temp_type);
  }

  fgets(the_line,sizeof(the_line), file1);

  pos = ftell(file1);

  while (fgets(the_line,sizeof(the_line), file1))
  {
    if (EQn(the_line,"FORMAT",6))
      break;

    if EQn(the_line,"CONECT",6)
    {
      the_line[80] = '\0';
      tokens = count_tokens(the_line," ");
      sscanf(&the_line[7],"%d",&i);
      for (j = 0; j < (tokens - 2); j++)
      {
	sscanf(&the_line[12 + j * 6],"%d",&con_atm);
	if  ((i <= Atoms) && (!bonded(mol,i,con_atm)))
	{
	  add_bond(mol,i,con_atm);
	}
      }
    }
  }

  fseek(file1,pos,0);
  for (i = 1; i <= Atoms; i++)
  {
    for (j = 0; j < Valence(i); j++)
    {
      BO(i,j) = 1;
    }
  }

  while (fgets(the_line,sizeof(the_line), file1))
  {
    if (EQn(the_line,"END",3))
      break;

    if EQn(the_line,"ORDER",5)
    {
      the_line[80] = '\0';
      tokens = count_tokens(the_line," ");
      sscanf(&the_line[7],"%d",&i);
      for (j = 0; j < (tokens - 2); j++)
      {
	sscanf(&the_line[12 + j * 6],"%d",&bond_order);
	BO(i,j) = bond_order;
      }
    }
  }
  
  build_connection_table(mol); 
  assign_types(mol);
  return(TRUE);
}








