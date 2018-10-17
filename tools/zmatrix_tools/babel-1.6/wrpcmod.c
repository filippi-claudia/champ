/****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact:

babel@mercury.aichem.arizona.edu
-------------------------------------------------------------------------------
FILE : wrpcmod.c
AUTHOR(S) : Abby Parrill
DATE : 6-94
PURPOSE : Write a pcmodel file into the UMS

*****/

#include "bbltyp.h"

int write_pcmod(FILE *file1, ums_type *mol)
{
  int i;
  int j;
  char temp_type[5];
  int type_name;
  int bo;
  
  fprintf(file1,"{PCM %s\n",OutfileName);
  fprintf(file1,"NA %d\n",Atoms);
  
  for (i = 1; i <= Atoms; i++)
  {
    get_output_type(i,"PCM",Type(i),temp_type,dummy);
    type_name = atoi(temp_type);
    type_name = update_mm2_types(mol,i,type_name);

    fprintf(file1,"AT %d %d %8.4f %8.4f %8.4f B", i, type_name, X(i), Y(i),
	    Z(i));
    for (j = 0; j < Valence(i); j++)
    {
      bo = get_bond_order(mol,i,Connection(i,j));
      fprintf(file1," %d %d", Connection(i,j), bo);
    }
    fprintf(file1,"\n");
  }
  fprintf(file1,"}\n");
  return(TRUE);
}


    

