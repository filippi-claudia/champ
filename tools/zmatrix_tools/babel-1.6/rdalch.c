/*****
  This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
  
  For more information please contact :
  
  babel@mercury.aichem.arizona.edu
  ----------------------------------------------------------------------------
  FILE : rdalch.c
  AUTHOR(S) : Pat Walters
  DATE : 12-92
  PURPOSE : routines to read an Alchemy file
  ******/

#include "bbltyp.h"

int read_alchemy(FILE *file1, ums_type *mol)
{
  int i;
  char input_line[BUFF_SIZE];
  char temp_type[5];
  char bo_string[10];
  int column;
  
  fgets(input_line,sizeof(input_line),file1);
  sscanf(input_line,"%d %*s %d",
	 &Atoms,
	 &Bonds);

  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);

  column = locate_input_type("ALC");
  for (i = 1; i <= Atoms; i ++)
  {
    UpdateProgress();
    fgets(input_line,sizeof(input_line),file1);
    sscanf(input_line,"%*d %s %lf %lf %lf",
	   temp_type,
	   &X(i),
	   &Y(i),
	   &Z(i));
    Atomic_number(i) = get_input_type(i,column,temp_type,Type(i),dummy);
  }
  for (i = 0; i < Bonds; i++)
  {
    fgets(input_line,sizeof(input_line),file1);
    sscanf(input_line,"%*d%d%d%s",&Start(i),&End(i),bo_string);
    Bond_order(i) = translate_alchemy_bond_order(bo_string);
  }
  dissect_connection_table(mol);
  return(TRUE);
}

int translate_alchemy_bond_order(char *bo_string)
{
  char err_string[50];
  
  if EQ(bo_string,"SINGLE")
    return(1);
  if EQ(bo_string,"DOUBLE")
    return(2);
  if EQ(bo_string,"TRIPLE")
    return(3);
  if EQ(bo_string,"AROMATIC")
    return(5);
  sprintf(err_string,"No bond type for Alchemy label %s",bo_string);
  return(1);
}



