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
PURPOSE : Routines to write format for use with conjure

******/
#include "bbltyp.h"


int 
write_conjure_tmplt(FILE *file1, ums_type *mol)
{ 
  int i,count;
  ums_type *new_ums = NULL;

  fprintf(file1,"%s\n",Title);
  fprintf(file1,"%9.4f\n",Energy);
  fprintf(file1,"%5d%5d\n",Atoms,Bonds);
  
  count = 1;
  for (i = 1;i <= Atoms;i++)
    if (NOTEQn("H",Type(i),1))
    {
      Redo(i) = count;
      count++;
    }
  
  for (i = 1;i <= Atoms;i++)
    if (EQn("H",Type(i),1))
    {
      Redo(i) = count;
      count++;
    }

  new_ums = (ums_type *)malloc(sizeof(ums_type));
  
  if (!new_ums)
  {
    show_warning("Unable to allocate memory for temporary ums");
    return(FALSE);
  }
  
  new_ums->num_atoms = Atoms;
  new_ums->num_bonds = Bonds;
  initialize_ums(&new_ums);
  strcpy(new_ums->title,Title);

  for (i = 1;i <= Atoms;i++)
  {
    new_ums->atoms[Redo(i)].point.x = X(i);
    new_ums->atoms[Redo(i)].point.y = Y(i);
    new_ums->atoms[Redo(i)].point.z = Z(i);
    strcpy(new_ums->atoms[Redo(i)].type,Type(i));
  }
  
  for (i = 0;i < Bonds;i++)
  {
    new_ums->connections[i].start = Redo(Start(i));
    new_ums->connections[i].end = Redo(End(i));
    new_ums->connections[i].bond_order = Redo(Bond_order(i));
  }
  
  dissect_connection_table(new_ums);

  mol = new_ums;
  
  for (i = 1;i <= Atoms;i++)
    qsort(mol->atoms[i].connected_atoms,Valence(i),
	  sizeof(int),QSORT_PROTO compare_int);
  
  
  for (i = 1;i <= Atoms; i++)
  {
    fprintf(file1,"%-3s %8.5f  %8.5f  %8.5f\n",
	    Type(i),
	    X(i),
	    Y(i),
	    Z(i));
  }

  build_connection_table(mol);

  for (i = 0;i < Bonds;i++)
    fprintf(file1,"%5d%5d%5d\n",Start(i),End(i),Bond_order(i));

  if (new_ums)
  {
    release_ums(new_ums);
    free(new_ums);
  }
  
  return(TRUE);
}











