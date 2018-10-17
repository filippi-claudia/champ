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

FILE : rdmacmod.c
AUTHOR(S) : Pat Walters
DATE : 3-93
PURPOSE :  Read a macromodel file

Modified 5-93 by Pat Walters
Allowed the use of single and multi-structure files

******/

#include "bbltyp.h"

int 
read_macromodel(FILE *file1, ums_type *mol)
{
  char mmd_line[BUFF_SIZE];
  int i;
  char temp_type[5];
  char chg_string[15];
  int column;
  

  fgets(mmd_line,sizeof(mmd_line),file1);
  sscanf(mmd_line,"%d",&Atoms);
  initialize_ums(&mol);
  sscanf(&mmd_line[47],"%lf",&Energy);
  Energy = Energy/4.184;
     
  ShowProgress(Atoms,"Reading Atoms");

  column = locate_input_type("MMD");
  for (i = 1; i <= Atoms; i++)
  {
    UpdateProgress();
    fgets(mmd_line,sizeof(mmd_line),file1);
    sscanf(mmd_line,"%s %d%d %d%d %d%d %d%d %d%d %d%d %lf%lf%lf",
	   temp_type,
	   &Connection(i,0),
	   &BO(i,0),
	   &Connection(i,1),
	   &BO(i,1),
	   &Connection(i,2),
	   &BO(i,2),
	   &Connection(i,3),
	   &BO(i,3),
	   &Connection(i,4),
	   &BO(i,5),
	   &Connection(i,5),
	   &BO(i,6),
	   &X(i),
	   &Y(i),
	   &Z(i));

    if (strlen(mmd_line) > 106)
    {
      my_strncpy(chg_string,&mmd_line[100],9);
      Charge(i) = atof(chg_string);
    }
    Atomic_number(i) = get_input_type(i,column,temp_type,Type(i),dummy);
  }  
  figure_valence(mol);
  build_connection_table(mol);

  return(TRUE);
}

void 
figure_valence(ums_type *mol)
{
  int i, j;
  
  for (i = MIN_ATOM; i <= Atoms; i++)
  {
    Valence(i) = 0;
    for (j = 0; j < 6; j++)
    {
      if (Connection(i,j) != 0)
	Valence(i)++;
    }
  }
}

int
assign_mmd_bond_order(ums_type *mol)
{
  int i,j,k;
  int found;
  
  for (i = 0; i < Bonds; i++)
  {
    found = FALSE;
    for (j = 0; j < Atoms; j++)
    {
      if (found == TRUE) break;
      for (k = 0; k < Valence(j); k++)
      {
	if ((Start(i) == j) && (End(i) == Connection(j,k)))
	{
	  Bond_order(i) = BO(j,k);
	  found = TRUE;
	  break;
	}
      }
    }
  }
  return(0);
}

    
    
	  


