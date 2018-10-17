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

FILE : renum.c
AUTHOR(S) : Pat Walters
DATE : 1-95
PURPOSE : routines to renumber a structure in order to produce a more reasonable
Z-matrix
	
******/

#include "bbltyp.h"

ums_type *renumber(ums_type *mol)
{
  vect_type v;
  
  center_at_origin(mol,&v);
  find_dist_from_origin(mol);
  sort_by_dist_to_origin(mol);
  mol = build_new_ums(mol,Atoms);
  return(mol);
}

void find_dist_from_origin(ums_type *mol)
{
  int i;
  coord_type origin;
  
  origin.x = 0.0;
  origin.y = 0.0;
  origin.z = 0.0;
  
  for (i = 1; i <= Atoms; i++)
  {
    Double(i) = distance(Point(i),origin);
  }
}

void sort_by_dist_to_origin(ums_type *mol)
{
  int i,j;
  temp_atom_rec *temp;

  temp = (temp_atom_rec *)malloc(Atoms * sizeof(temp_atom_rec));

  if (!temp)
    fatal_error("Error allocating memory in sort_by_dist_to_origin");

  for (i = 0; i < Atoms; i++)
  {
    j = i + 1;
    temp[i].x = X(j);
    temp[i].y = Y(j);
    temp[i].z = Z(j);
    temp[i].num = j;
    temp[i].dist = Double(j);
  }
  
  qsort(temp,Atoms,sizeof(temp_atom_rec),QSORT_PROTO sort_by_dist);

  for (i = 0; i < Atoms; i++)
  {
    printf("%d %10.3f%10.3f%10.3f - %10.3f\n",
	   temp[i].num,temp[i].x,temp[i].y,temp[i].z,temp[i].dist);
  }
  
  j = 1;
  for (i = 0; i < Atoms; i++)
  {
    if (Type(temp[i].num)[0] != 'H')
    {
      Redo(temp[i].num) = j;
      j++;
    }
  }
  for (i = 0; i < Atoms; i++)
  {
    if (Type(temp[i].num)[0] == 'H')
    {
      Redo(temp[i].num) = j;
      j++;
    }
  }
}


int sort_by_dist(temp_atom_rec *a, temp_atom_rec *b)
{
  if (a->dist > b->dist)
    return(1);
  if (a->dist < b->dist)
    return(-1);
  return(0);
}


