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

FILE : rdboogie
AUTHOR(S) : Matt Stahl
DATE : 7-93
PURPOSE : Routines to read a DISCO output file

******/

#include "bbltyp.h"

int 
read_boogie(FILE *file1, ums_type *mol)
{
  char boog_line[BUFF_SIZE];
  int i = 0,j;
  int count = 0;
  int pos_found = FALSE;
  long pos = 0;
  int result;
    
  while (fgets(boog_line,sizeof(boog_line), file1) != NULL)
  {
    if (memcmp(boog_line,"*",1) == 0)
    {
      if (pos_found == FALSE)
      {
	pos_found = TRUE;
	pos = ftell(file1);
      }
      count ++;
    }
  }

  Atoms = count;

  ShowProgress(Atoms,"Reading Atoms");

  initialize_ums(&mol);
  
  fseek(file1,pos,0);

  i = 1;
  
  while (fgets(boog_line,sizeof(boog_line), file1) != NULL)
  {
    if (memcmp(boog_line,"*",1) != 0)
    {
      UpdateProgress();
	j = sscanf(boog_line,"%s %lf %lf %lf",Type(i),&X(i),&Y(i),&Z(i)); 
	if (j != 4)
	{
	  show_warning("Input file error");
	  release_ums(mol);
	  Atoms = Bonds = 0;
	  return(FALSE);
	}
	clean_atom_type(Type(i));
	i ++;
      }
  }
  result = assign_radii(mol);
  result = assign_bonds(mol);
  result = assign_types(mol);
  result = build_connection_table(mol);
  assign_bond_order(mol);

  return(TRUE);
}

       
























