/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : precip.c
AUTHOR(S) : Pat Walters
DATE : 8-23-95
PURPOSE : remove atoms with valence 0 from a UMS

******/

#include "bbltyp.h"

ums_type *precipitate(ums_type *mol)
 {  
  ums_type *head;

  head = dissect_ums(mol);
  cleanup_ums_lst(head->next);
  return(head);
}

int tag_salts(ums_type *mol)
{
  int i;
  int j = 0;

  for (i = 1;i <= Atoms;i++)
  {   
    if (Valence(i) == 0)
    {
      Redo(i) = 0;
    }
    else
    {
      j++;
      Redo(i) = j;
    }
  }
  return(j);
}

void copy_ums(ums_type *dest, ums_type *src)
{
  int i,j;
  
  dest->num_atoms = src->num_atoms;
  dest->num_bonds = src->num_bonds;
  initialize_ums(&dest);
  strcpy(dest->title,src->title);
  
  for (i = 1; i <= src->num_atoms; i++)
  {
    strcpy(dest->atoms[i].type,src->atoms[i].type);
    dest->atoms[i].point.x = src->atoms[i].point.x;
    dest->atoms[i].point.y = src->atoms[i].point.y;
    dest->atoms[i].point.z = src->atoms[i].point.z;
    dest->atoms[i].max_bonds = src->atoms[i].max_bonds;
    dest->atoms[i].valence = src->atoms[i].valence;
    dest->atoms[i].redo = 0;
    for (j = 0; j < src->atoms[i].valence; j++)
    {
      dest->atoms[i].connected_atoms[j] = src->atoms[i].connected_atoms[j];
      dest->atoms[i].bond_order[j] = src->atoms[i].bond_order[j];
    }
    dest->atoms[i].radius = src->atoms[i].radius;
    dest->atoms[i].bond_ord_rad = src->atoms[i].bond_ord_rad;
    dest->atoms[i].dble = src->atoms[i].dble;
    dest->atoms[i].organic = src->atoms[i].organic;
    dest->atoms[i].charge = src->atoms[i].charge;
    dest->atoms[i].redo = src->atoms[i].redo;
  }
  for (i = 0; i < src->num_bonds; i++)
  {
    dest->connections[i].start = src->connections[i].start;
    dest->connections[i].end = src->connections[i].end;
    dest->connections[i].bond_order = src->connections[i].bond_order;
  }
}

    







