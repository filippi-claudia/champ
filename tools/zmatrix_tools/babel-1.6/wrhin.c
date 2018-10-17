/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : wrhin.c
AUTHOR(S) : Pat Walters
DATE : 1-94
PURPOSE : Routines to write a Hyperchem HIN file

******/
#include "bbltyp.h"

static int file_num = 1;

int 
write_hin(FILE *file1, ums_type *mol)
{ 
  int i,j;
  char type_name[5];
  int result;
  char ord_sym;
  int bo;

  fprintf(file1,"mol %d %s\n",file_num,OutfileName);
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1,"atom %d - %-3s **  - %8.5f %8.5f  %8.5f  %8.5f %d ",
	    i,
	    type_name,
	    Charge(i),
	    X(i),
	    Y(i),
	    Z(i),
	    Valence(i));
    for (j = 0; j < Valence(i); j++)
    {
      bo = get_bond_order(mol,i,Connection(i,j));
      switch(bo)
      {
      case 2 :
	ord_sym = 'd';
	break;
      case 3 :
	ord_sym = 't';
	break;
      default :
	ord_sym = 's';
	break;
      }
      fprintf(file1,"%d %c ",Connection(i,j),ord_sym);
    }
    fprintf(file1,"\n");
  }
  fprintf(file1,"endmol %d\n",file_num);
  file_num++;
  return(TRUE);
}


  
      









