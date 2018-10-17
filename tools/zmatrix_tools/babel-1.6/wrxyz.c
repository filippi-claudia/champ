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

FILE : wrxyz.c
AUTHOR(S) : Pat Walters
DATE : 11-92
PURPOSE : Routines to write an XYZ file as used by the Xmol program from MSC
******/

#include "bbltyp.h"

int 
write_xyz(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  int result;

  fprintf(file1,"%d\n",Atoms);
  fprintf(file1,"%s\t%15.7f\n",Title,Energy);
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1,"%3s%15.5f%15.5f%15.5f\n",
	    type_name,
	    X(i),
	    Y(i),
	    Z(i));
  }
  return(TRUE);
}


  
      









