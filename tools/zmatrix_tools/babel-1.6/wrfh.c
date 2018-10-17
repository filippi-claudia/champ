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

FILE : wrfh.c
AUTHOR(S) : Matthew Stahl
DATE : 5-94
PURPOSE : Routines to write a fenske hall internal coordinate file, also contains 
internal to cartesian coordinate conversion routines.  This module was stolen directly
from some of Pat's mopac internal coordinate file routines.  :) 

******/

#include "bbltyp.h"

int 
  write_fenske_zmat(FILE *file1, ums_type *mol)
{
  int i=0;
  char type_name[5];
  int result;
  
  
  if (mol->internal == NULL)
  {
    initialize_internal(&mol);
    cartint(mol);
    cartgeom(mol);
  }

  fprintf(file1,"%s\n",OutputKeywords);
  fprintf(file1,"%d\n",Atoms);
    
  for (i = 1;i <= Atoms;i++)
  { 
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    if (i == 1)
    {
      fprintf(file1,"%-2s 1\n",type_name);
    }
    if (i == 2)
    {
      fprintf(file1,"%-2s%3d%6.3f\n",
	      type_name,
	      mol->internal[i].na,
	      mol->internal[i].r); 
    }
    if (i == 3)
    {
      fprintf(file1,"%-2s%3d%6.3f%3d%8.3f\n",
	      type_name,
	      mol->internal[i].na,
	      mol->internal[i].r,
	      mol->internal[i].nb,
	      mol->internal[i].w); 

    }
    if (i > 3)
    {
      if (mol->internal[i].t < 0)
	mol->internal[i].t += 360;
      
      fprintf(file1,"%-2s%3d%6.3f%3d%8.3f%3d%6.1f\n",
	      type_name,
	      mol->internal[i].na,
	      mol->internal[i].r,
	      mol->internal[i].nb,
	      mol->internal[i].w,
	      mol->internal[i].nc,
	      mol->internal[i].t);
    }
  }
  return(TRUE);
}













