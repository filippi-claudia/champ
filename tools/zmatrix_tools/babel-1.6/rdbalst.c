/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
-----------------------------------------------------------------------------
FILE : rdbalst.c
AUTHOR(S) : Pat Walters
DATE : 10-16-93
PURPOSE : routines to read Ball and Stick format files
******/

#include "bbltyp.h"

int 
read_bs(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i,j;
  int result;
  int tokens;
  
  fgets(the_line,sizeof(the_line), file1);
  fgets(the_line,sizeof(the_line), file1);
  sscanf(the_line,"%d",&Atoms);

  ShowProgress(Atoms,"Reading Atoms");

  result = initialize_ums(&mol);
  for (i = MIN_ATOM; i <= Atoms;i ++)
    {

      UpdateProgress();

      fgets(the_line,sizeof(the_line), file1);
      tokens = count_tokens(the_line,"\t\n ");
      Valence(i) = tokens - 4;
      sscanf(the_line,"%s %lf %lf %lf",
	     Type(i),
	     &X(i),
	     &Y(i),
	     &Z(i));
      for (j = 0; j < Valence(i); j++)
      {
	Connection(i,j) = atoi(gettoken(the_line,"\t\n ",j + 5));
      }
      clean_atom_type(Type(i));
    }

  result = assign_radii(mol);
  result = assign_types(mol);
  result = build_connection_table(mol);
  assign_bond_order(mol);
  return(TRUE);
}







   
    
    
    
	  


