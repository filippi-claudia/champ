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

FILE : wrmicro.c
AUTHOR(S) : Pat Walters
DATE : 10-27-93
PURPOSE : Routines to write a Micro World file

******/



#include "bbltyp.h"

int 
write_micro(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  
  fprintf(file1,"\n");
  
  for(i = 1;i <= Atoms; i++)
  {
    get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1,"%5d%-3s %8.5f  %8.5f  %8.5f \n",
	    i,
	    type_name,
	    X(i),
	    Y(i),
	    Z(i));
  }
  return(TRUE);
}


  
      









