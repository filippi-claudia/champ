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

FILE : rdmicro.c
AUTHOR(S) : Pat Walters
DATE : 10-93
PURPOSE : routines to read a MicroWorld
******/

#include "bbltyp.h"


int 
read_micro(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i = 0;
  int result;

  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (count_tokens(the_line,"\n\t ") >= 4)
      i++;
  }
  Atoms = i;
  rewind(file1);
#ifdef MAC
ShowProgress(Atoms,"Reading Atoms");
#endif
  result = initialize_ums(&mol);
  i = 1;
  while ( i <= Atoms)
  {
    fgets(the_line,sizeof(the_line), file1);
    if (count_tokens(the_line,"\t\n ") == 4)
    {
#ifdef MAC
      UpdateProgress();
#endif
      sscanf(the_line,"%s%lf%lf%lf",Type(i),&X(i),&Y(i),&Z(i));
      strcpy(Type(i), strip_front_num(Type(i)));
      i ++;
    }
  }
  if (Atoms == 0)
  {
    return(FALSE);
  }
  result = assign_radii(mol);
  result = assign_bonds(mol);
  result = assign_types(mol);
  result = build_connection_table(mol);
  assign_bond_order(mol);
  return(TRUE);  
}



char *strip_front_num(char *the_str)
{
  while (isdigit(the_str[0]))
    the_str++;
  return(the_str);
}
    
	  


