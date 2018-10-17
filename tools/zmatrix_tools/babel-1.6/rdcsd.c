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

FILE : rdcsd.c
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : routines to read a CSD GSTAT file

******/

#include "bbltyp.h"

int 
read_csd(FILE *file1, ums_type *mol)
{
  char csd_line[BUFF_SIZE];
  int i = 0;
  int result;
  long pos;
  
  fgets(csd_line,sizeof(csd_line), file1);
  if (strstr(csd_line,"FRAG") != NULL)
    {
      sscanf(csd_line,"%s",Title);
      Title[8] = '\0';
      pos = ftell(file1);
      i = 0;

      while ((fgets(csd_line,sizeof(csd_line), file1) != NULL) 
	     && (strstr(csd_line,"FRAG") == NULL)) 
	{
	  if (count_tokens(csd_line," ") == 4)
	    i++;
	}
      Atoms = i;
      fseek(file1,pos,0);
ShowProgress(Atoms,"Reading Atoms");
      result = initialize_ums(&mol);
      
      i = 1;
      while (i <= Atoms)
	{
	  UpdateProgress();
	  fgets(csd_line,sizeof(csd_line), file1);
	  if (count_tokens(csd_line," ") == 4)
	    {
	      sscanf(csd_line,"%s%lf%lf%lf",Type(i),&X(i),&Y(i),&Z(i));
	      clean_atom_type(Type(i));
	      i ++;
	    }
	}
	result = assign_radii(mol);
	result = assign_bonds(mol);
	result = assign_types(mol);
	result = build_connection_table(mol);
      assign_bond_order(mol);
  }
  else
    return(FALSE);
  return(TRUE);  
}

   
    
    
    
	  


