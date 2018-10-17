/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : wrbalst.c
AUTHOR(S) : Pat Walters
DATE : 10-10-93
PURPOSE : Routines to write a Tinker XYZ file

******/
#include "bbltyp.h"

int 
write_tinker(FILE *file1, ums_type *mol)
{ 
  int i,j;
  char xyz_name[10];
  char mm2_name[10];
  int result;
  int type_name;
  
  fprintf(file1,"%6d %-20s\n",Atoms,Title);
  
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),xyz_name,all_caps);
    result = get_output_type(i,"MM2",Type(i),mm2_name,dummy);
    type_name = atoi(mm2_name);
    type_name = update_mm2_types(mol,i,type_name);
    fprintf(file1,"%6d %2s  %12.6f%12.6f%12.6f %5d",
	    i,
	    xyz_name,
	    X(i),
	    Y(i),
	    Z(i),
	    type_name);
    for (j = 0; j < Valence(i); j++)
      fprintf(file1,"%6d",Connection(i,j));
    fprintf(file1,"\n");
  }
  return(TRUE);
}


  
      









