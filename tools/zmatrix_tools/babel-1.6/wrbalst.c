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

FILE : wrbalst.c
AUTHOR(S) : Pat Walters
DATE : 10-10-93
PURPOSE : Routines to write a Ball and Stick file

******/
#include "bbltyp.h"

int 
write_bs(FILE *file1, ums_type *mol)
{ 
  int i,j;
  char type_name[5];
  int result;
  
  fprintf(file1,"%s\n",Title);
  fprintf(file1,"%d\n",Atoms);
  
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1,"%-3s %8.5f  %8.5f  %8.5f",
	    type_name,
	    X(i),
	    Y(i),
	    Z(i));
    for (j = 0; j < Valence(i); j++)
      fprintf(file1,"%6d",Connection(i,j));
    fprintf(file1,"\n");
  }
  return(TRUE);
}


  
      









