/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
------------------------------------------------------------------------------

FILE : wrunichm.c
AUTHOR(S) : Pat Walters
DATE : 2-95
PURPOSE : Routines to write an Cray UniChem file
******/

#include "bbltyp.h"

int 
write_unichem(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  int result;

  fprintf(file1,"%s\n",Title);
  fprintf(file1,"%d\n",Atoms);
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1,"%3d%15.5f%15.5f%15.5f\n",
	    Atomic_number(i),
	    X(i),
	    Y(i),
	    Z(i));
  }
  return(TRUE);
}

  
      









