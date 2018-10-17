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

FILE : wrwiz.c
AUTHOR(S) : Matthew Stahl
DATE : 5-94
PURPOSE : Routines to write format for use with Wizard

******/
#include "bbltyp.h"


int 
write_wizard(FILE *file1, ums_type *mol)
{ 
  int i;
  ums_type *new_ums = NULL;
  
  fprintf(file1,"%s\n",Title);
  fprintf(file1,"%9.4f\n",Energy);
  fprintf(file1,"%5d%5d\n",Atoms,Bonds);
  
  for (i = 1;i <= Atoms;i++)
    qsort(mol->atoms[i].connected_atoms,Valence(i),
	  sizeof(int),QSORT_PROTO compare_int);

  for (i = 1;i <= Atoms; i++)
  {
    if (EQ(Type(i),"O-"))
    {
      if (Valence(i) == 1)
	if (EQ(Type(Connection(i,0)),"Cac"))
	  strcpy(Type(i),"O2");

      if (Valence(i) == 2)
	if (EQ(Type(Connection(i,0)),"Cac") || 
	    EQ(Type(Connection(i,0)),"Cac"))
	  strcpy(Type(i),"O3");
    }
  }
  
  for (i = 1;i <= Atoms; i++)
  {
    fprintf(file1,"%-3s %8.5f  %8.5f  %8.5f\n",
	    Type(i),
	    X(i),
	    Y(i),
	    Z(i));
  }

  for (i = 0;i < Bonds;i++)
    fprintf(file1,"%5d%5d%5d\n",Start(i),End(i),Bond_order(i));

  if (new_ums)
    release_ums(new_ums);
  
  return(TRUE);
}


/*
int compare_int(int *a,int *b)
{
  int result=0;
  
  if (*a < *b) result = -1;
  else if
    (*a == *b) result = 0;
  else if
    (*a > *b) result =  1;
  
  return(result);
}
*/

int compare_int(int *a, int *b)
{
  return(*a < *b ? -1 : *a != *b);
}











