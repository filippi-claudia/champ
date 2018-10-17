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
FILE : rdsybmol.c
AUTHOR(S) : Pat Walters
DATE : 2-94
PURPOSE : routines to a sybyl MOL file
******/

#include "bbltyp.h"

int 
read_sybyl_mol(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i;
  int result;
  char temp1[5],temp2[15],temp3[5];
  int column;
  
  fgets(the_line,sizeof(the_line), file1);
  sscanf(the_line,"%d",&Atoms);
  ShowProgress(Atoms,"Reading Atoms");
  result = initialize_ums(&mol);

  column = locate_input_type("MOL");
  for (i = MIN_ATOM; i <= Atoms;i++)
  {
    UpdateProgress();
    fgets(the_line,sizeof(the_line), file1);
    sscanf(the_line,"%*d %s %lf %lf",
	   temp1,
	   &X(i),
	   &Y(i));
    my_strncpy(temp2,&the_line[26],9);
    my_strncpy(temp3,&the_line[35],3);
    Z(i) = atof(temp2);
    Atomic_number(i) = get_input_type(i,column,temp1,Type(i),dummy);    
  }
  fgets(the_line,sizeof(the_line), file1);
  sscanf(the_line,"%d",&Bonds);
  for (i = 0; i < Bonds; i++)
  {
    fgets(the_line,sizeof(the_line), file1);
    sscanf(the_line,"%*d %d %d %d",
	   &Start(i),
	   &End(i),
	   &Bond_order(i));
  }
  dissect_connection_table(mol);
  fgets(the_line,sizeof(the_line), file1);
  return(TRUE);
}
