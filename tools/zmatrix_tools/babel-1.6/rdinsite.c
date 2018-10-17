/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
-------------------------------------------------------------------------------
FILE : rdbalst.c
AUTHOR(S) : Pat Walters
DATE : 10-16-93
PURPOSE : routines to read Ball and Stick format files
******/

#include "bbltyp.h"

int 
read_biosym_car(FILE *file1, ums_type *mol)
{
  char the_line[300];
  char pbc_line[300];
  int i;
  int result;
  long pos;
  int done = FALSE;
  int count = 3;
  
  while ((fgets(the_line,sizeof(the_line),file1) != NULL) && (!done))
  {
  	if (strstr(the_line,"PBC"))
	{
	  done = TRUE;
	  strcpy(pbc_line,the_line);
	}
  }
  if (strstr(pbc_line,"ON"))
  	count = 2;
  else
  	count = 1;
  for (i = 0; i < count; i++)
  	fgets(the_line,sizeof(the_line),file1);
  pos = ftell(file1);
  Atoms = -1;
  the_line[0] = '\0';
  while (strstr(the_line,"end") == NULL)
  {
    fgets(the_line,sizeof(the_line),file1);
    Atoms++;
  }
  fseek(file1,pos,0);
#ifdef MAC
  ShowProgress(Atoms,"Reading Atoms");
#endif
  result = initialize_ums(&mol);
  for (i = MIN_ATOM; i <= Atoms;i ++)
  {
#ifdef MAC
    UpdateProgress();
#endif
    fgets(the_line,sizeof(the_line), file1);
    sscanf(the_line,"%*s %lf %lf %lf %*s%*s%*s %s",
	   &X(i),
	   &Y(i),
	   &Z(i),
	   Type(i));
    clean_atom_type(Type(i));
  }
  fgets(the_line,sizeof(the_line),file1);
  fgets(the_line,sizeof(the_line),file1);
  result = assign_radii(mol);
  result = assign_bonds(mol);
  result = assign_types(mol);
  result = build_connection_table(mol);
  assign_bond_order(mol);
  return(TRUE);
}







   
    
    
    
	  


