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

FILE : rdcharmm.c
AUTHOR(S) : Pat Walters
DATE : 12-93
PURPOSE : read a CHARMm .CRD file

******/
#include "bbltyp.h"

int 
  read_charmm(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int result;
  int done = FALSE;
  int i;
  
  while ((fgets(the_line,sizeof(the_line), file1) != NULL) && (the_line[0] == '*'));
  sscanf(the_line,"%d",&Atoms);
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  initialize_residues(&mol);

  for (i = 1; i <= Atoms; i++)
  {
    fgets(the_line,sizeof(the_line), file1);
    UpdateProgress();
    ChainNum(i) = 0;
    sscanf(the_line,"%*s %d%s%s",&ResNum(i),ResName(i),AtmId(i));

    if (isspace(the_line[15]))
    {
      Type(i)[0] = the_line[16];
      Type(i)[1] = '\0';
    }
    else
    {
      Type(i)[0] = the_line[15];
      Type(i)[1] = the_line[16];
      Type(i)[2] = '\0';
    }
    sscanf(&the_line[22],"%lf",&X(i));
    sscanf(&the_line[32],"%lf",&Y(i));
    sscanf(&the_line[42],"%lf",&Z(i));
  }
  result = assign_radii(mol);  
  result = assign_bonds(mol);
  result = assign_types(mol);
  result = build_connection_table(mol);
  assign_bond_order(mol);
  return(TRUE);  
}

   
    
    
    
	  


