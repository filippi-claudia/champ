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
DATE : 12-94
PURPOSE : Routines to write an icon8 input file
******/

#include "bbltyp.h"

int 
write_icon8(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  int result;

  fprintf(file1,"%s\n",Title);
  fprintf(file1,"%6d%12s",Atoms,"");
  fprintf(file1,"FFFFF 0.00 0.000 0.000FTTTTTTFFTTTTTTFTTTTFFFFFFFFFTTFFFFFFFFF\n");
  for (i = 1;i <= Atoms; i++)
  { 
     fprintf(file1,"%15.6f%15.6f%15.6f \n",
	     X(i),
	     Y(i),
	     Z(i));
  }

  for (i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1,"%2s",type_name);
    if (((i % 40) == 0 ) && (i != 0))
      fprintf(file1,"\n");
  }
  fprintf(file1,"\n");
  return(TRUE);
}


  
      









