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

FILE : wrchdrw.c
AUTHOR(S) : Pat Walters
DATE : 1-10-92
PURPOSE : Routines to write a ChemDraw connection table
******/

#include "bbltyp.h"

int 
write_chem_draw(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  int result;

  char axis;
/*  axis = OutputKeywords[0]; */
  axis = 'x';
  
  
  fprintf(file1,"%s\n",Title);
  fprintf(file1," %d %d\n",Atoms,Bonds);
  
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    switch(axis)
    {
    case 'x':
    fprintf(file1," %9.4f %9.4f    0.0000 %-1s\n",
	    Y(i),
	    Z(i),
	    type_name);
    break;
    case 'y':
    fprintf(file1," %9.4f %9.4f    0.0000 %-1s\n",
	    X(i),
	    Z(i),
	    type_name);
    break;
    case 'z':
    fprintf(file1," %9.4f %9.4f    0.0000 %-1s\n",
	    X(i),
	    Y(i),
	    type_name);
    break;
    default :
    {
    fprintf(file1," %9.4f %9.4f %9.4f %-5s\n",
	    X(i),
	    Y(i),
	    Z(i),
	    type_name);
    }
  }
  }

  for(i = 0;i < Bonds; i++)
  {
    fprintf(file1,"%3d%3d%3d%3d\n",
	    Start(i),
	    End(i),
	    Bond_order(i),1);
  }
  return(TRUE);
}















