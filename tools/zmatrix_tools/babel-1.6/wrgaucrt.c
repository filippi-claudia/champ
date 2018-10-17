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

FILE : wrxyz.c
AUTHOR(S) : Pat Walters
DATE : 1-94
PURPOSE : Routines to write a Gaussian Cartesian file
******/

#include "bbltyp.h"

int 
write_gaus_crt(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  int result;
  
  fprintf(file1,"%cmem=20000000\n",'\045');
  fprintf(file1,"#%s\n\n",OutputKeywords);
  
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    clean_atom_type(type_name);
    fprintf(file1,"%-3s0      x%-5d     y%-5d     z%-5d \n",
	    type_name,i,i,i);
  }
  for (i = 1; i <= Atoms; i++)
    fprintf(file1,"x%-4d %10.5f\n",i,X(i));
  for (i = 1; i <= Atoms; i++)
    fprintf(file1,"y%-4d %10.5f\n",i,Y(i));
  for (i = 1; i <= Atoms; i++)
    fprintf(file1,"z%-4d %10.5f\n",i,Z(i));
  return(TRUE);
}


  
      









