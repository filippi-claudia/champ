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

FILE : wrsybmol.c
AUTHOR(S) : Pat Walters
DATE : 2-94
PURPOSE : Routines to write a Sybyl MOL file
******/

#include "bbltyp.h"

int 
write_sybyl_mol(FILE *file1, ums_type *mol)
{ 
  int i;
  int type_name;
  char temp_type[5];
  char title_str[20];
  char ele[3];

  strncpy(title_str,OutfileName,20);
  
  fprintf(file1,"%4d MOL",Atoms);
  fprintf(file1,"%20s%11s%d\n",OutfileName,"",0);

  for(i = 1;i <= Atoms; i++)
  {
    get_element_type(mol,i,ele);
    get_output_type(i,"MOL",Type(i),temp_type,dummy);
    type_name = atoi(temp_type);

    fprintf(file1,"%4d%4d%9.4f%9.4f%9.4f%s\n",
	    i,
	    type_name,
	    X(i),
	    Y(i),
	    Z(i),
	    ele);
  }
  fprintf(file1,"%4d MOL\n",Bonds);
  for(i = 0;i < Bonds; i++)
  {
    fprintf(file1,"%4d%4d%4d%9s%4d\n",
	    i + 1,
	    Start(i),
	    End(i),
	    "",
	    Bond_order(i));
  }
  fprintf(file1,"%4d MOL\n",0);
  return(TRUE);
}

