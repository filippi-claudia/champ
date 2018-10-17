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

FILE : wrcacao
AUTHOR(S) : Matt Stahl
DATE : 7-93
PURPOSE : Routines to write an cartesian coordinate file as used by the
          program cacao

******/

#include "bbltyp.h"

int
write_caccrt(FILE *file1, ums_type *mol)
{
  int i;
  char type_name[5];
  int result;

  fprintf(file1,"%s\n",Title);
  fprintf(file1,"%3d   DIST  0  0  0\n",Atoms);
  fprintf(file1,"CELL 1.,1.,1.,90.,90.,90.\n");

  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,dummy);
    fprintf(file1,"%2s %7.4f, %7.4f, %7.4f \n",
            type_name,
            X(i),
            Y(i),
            Z(i));
  }
  return(TRUE);
}



