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

FILE : wralch.c
AUTHOR(S) : Pat Walters
DATE : 12-92
PURPOSE : Write an alchemy file

******/

#include "bbltyp.h"

int
  write_alchemy(FILE *file1, ums_type *mol)
{
  int i;
  char temp_type[5];
  int result;
  char bond_string[10];
  
  fprintf(file1,"%5d ATOMS, %5d BONDS,     0 CHARGES\n",
	  Atoms,
	  Bonds);
  
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"ALC",Type(i),temp_type,dummy);
    fprintf(file1,"%5d %-6s%8.4f %8.4f %8.4f     0.0000\n",
	    i,
	    temp_type,
	    X(i),
	    Y(i),
	    Z(i));
  }

  for(i = 0; i < Bonds; i++)
  {
    switch(Bond_order(i))
    {
    case 1 :
      strcpy(bond_string,"SINGLE");
      break;
    case 2 :
      strcpy(bond_string,"DOUBLE");
      break;
    case 3 :
      strcpy(bond_string,"TRIPLE");
      break;
    case 5 :
      strcpy(bond_string,"AROMATIC");
      break;
    default :
      strcpy(bond_string,"SINGLE");
    }
    fprintf(file1,"%5d  %4d  %4d  %s\n",
	    i + 1,
	    Start(i),
	    End(i),
	    bond_string);
  }
  return(1);
}




void strip_extension(char *file_name, char *new_name)
{
  char  *pos;
  int end_pos;
  
  pos = strchr(file_name,'.');
  if (pos != NULL)
  {
    end_pos = pos - file_name;
    strncpy(new_name,file_name,pos - file_name);
    new_name[end_pos] = '\0';
  }
  else
    strcpy(new_name,file_name);
}






