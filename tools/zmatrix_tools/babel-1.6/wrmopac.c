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

FILE : wrmopac.c
AUTHOR(S) : Pat Walters (slightly modified by A. Parrill 7/94)
DATE : 11-92
PURPOSE : Routines to write a MOPAC cartesian coordinate input file

******/

#include "bbltyp.h"

int 
write_mopac(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  
  fprintf(file1,"%s\n",OutputKeywords);
  fprintf(file1,"%s\n",OutfileName);
  fprintf(file1,"%s\n",Title);
  
  for(i = 1;i <= Atoms; i++)
  {
    get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1,"%-3s%8.5f 1 %8.5f 1 %8.5f 1\n",
	    type_name,
	    X(i),
	    Y(i),
	    Z(i));
  }
  return(TRUE);
}

  
  
      









