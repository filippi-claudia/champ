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

FILE : wrmdl.c
AUTHOR(S) : Pat Walters
DATE : 12-01-93
PURPOSE : Routines to write an MDL molfile
******/

#include "bbltyp.h"

int 
write_molfile(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  int result;
  
  fprintf(file1,"%s\n",Title);
  fprintf(file1,"%s\n\n","  -ISIS-            3D");
  fprintf(file1,"%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d\n",
	  Atoms,Bonds,0,0,0,0,0,0,0,0,0);
  
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1,"%10.4f%10.4f%10.4f %-3s%2d%3d%3d%3d%3d\n",
	    X(i),
	    Y(i),
	    Z(i),
	    type_name,
	    0,0,0,0,0);
  }
  for(i = 0; i < Bonds; i++)
  {
    if (Bond_order(i) == 5)
      Bond_order(i) = 4;
    fprintf(file1,"%3d%3d%3d%3d%3d%3d\n",
	    Start(i),
	    End(i),
	    Bond_order(i),
	    0,0,0);
  }
  fprintf(file1,"M  END\n");
  return(TRUE);
}


  
      










