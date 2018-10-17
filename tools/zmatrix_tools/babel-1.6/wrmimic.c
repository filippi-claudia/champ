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

FILE : wrmimic.c
AUTHOR(S) : Pat Walters
DATE : 5-10-93
PURPOSE : Routines to write a MacMimic file

******/

#include "bbltyp.h"

int 
write_mac_mimic(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  
  fprintf(file1,"*Molecule\n");
  fprintf(file1,"%5d%5d \n",Atoms,Bonds);
  fprintf(file1,"*Model\nS\n*Display\nFFF\n");
  fprintf(file1,"*Position\n 200  200\n");

  fprintf(file1,"*Atom\n");  
  for(i = 1;i <= Atoms; i++)
  {
    get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1,"%-4d  %-3s %8.4f  %8.4f  %8.4f %d 0  \n",
	    i,
	    type_name,
	    X(i),
	    Y(i),
	    Z(i),
	    Valence(i));
  }

  fprintf(file1,"*Bond\n");  
  for(i = 0;i < Bonds; i++)
  {
    fprintf(file1,"%5d%5d\n",
	    Start(i),
	    End(i));
  }
  fprintf(file1,"*End\n");  
  return(TRUE);
}

