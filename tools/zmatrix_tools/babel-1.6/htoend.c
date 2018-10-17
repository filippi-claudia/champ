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

FILE : htoend.c
AUTHOR(S) : Matthew Stahl
DATE : 10-95
PURPOSE : Shift hydrogens to the end of the file

******/

#include "bbltyp.h"

void shift_h_to_end(ums_type *mol)
{
  int i,count;
  ums_type *new_ums = NULL;

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
    fatal_error("Unable to allocate memory for temporary ums");
  
  new_ums->num_atoms = Atoms;
  new_ums->num_bonds = Bonds;
  initialize_ums(&new_ums);
  
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
  
  release_ums(mol);
  mol = new_ums;
  
  for (i = 1;i <= Atoms;i++)
    qsort(mol->atoms[i].connected_atoms,Valence(i),
	  sizeof(int),QSORT_PROTO compare_int);
  
/*  if (new_ums)
    release_ums(new_ums);*/
}


