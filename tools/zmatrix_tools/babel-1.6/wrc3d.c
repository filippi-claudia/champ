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

FILE : wrc3d.c
AUTHOR(S) : Pat Walters
DATE : 10-19-93
PURPOSE : Routines to write a Chem 3d file

******/
#include "bbltyp.h"

int 
write_chem3d2(FILE *file1, ums_type *mol)
{
  write_chem3d(file1,mol,"C3D");
  return(TRUE);
}

int 
write_chem3d1(FILE *file1, ums_type *mol)
{
  write_chem3d(file1,mol,"MM2");
  return(TRUE);
}

int 
write_mmads(FILE *file1, ums_type *mol)
{
  write_chem3d(file1,mol,"MMADS");
  return(TRUE);
}

int 
write_chem3d(FILE *file1, ums_type *mol, char *mol_typ)
{ 
  int i,j;
  char type_name[5];
  int result;
  char ele_type[5];
  int atnum;
  int type_num;

  fprintf(file1,"%d",Atoms);
  if (EQ(mol_typ,"MMADS"))
  {
    fprintf(file1," %s",Title);
    strcpy(mol_typ,"MM2");
  }
  fprintf(file1,"\n");

  for(i = 1;i <= Atoms; i++)
  {
    result = xlate_std_type(mol_typ,Type(i),type_name);
    if (result == 0)
    {
      fprintf(stderr,"Unable to assign %s type to atom %d type = %s\n",
	      mol_typ,i,Type(i));
      atnum = get_atomic_number(Type(i));
      type_num = atnum * 10 + Valence(i);
      sprintf(type_name,"%d",type_num);
    }
    get_element_type(mol,i,ele_type);
    fprintf(file1,"%-3s %-5d %8.5f  %8.5f  %8.5f %5s",
	    ele_type,
	    i,
	    X(i),
	    Y(i),
	    Z(i),
	    type_name);
    for (j = 0; j < Valence(i); j++)
      fprintf(file1,"%6d",Connection(i,j));
    fprintf(file1,"\n");
  }
  return(TRUE);
}


      









