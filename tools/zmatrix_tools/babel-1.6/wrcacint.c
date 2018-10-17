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

FILE : wrcacint.c
AUTHOR(S) : Pat Walters
DATE : 7-94
PURPOSE : Routines to write a CACAO internal coordinate file

******/

#include "bbltyp.h"

int write_cacao_internal(FILE *file1, ums_type *mol)
{
  int i;
  char the_type[5];

  add_dummy_atoms(mol);
  center_at_first_atom(mol);
  initialize_internal(&mol);
  set_hilderbrandt_connections(mol);
  set_hilderbrandt_geometry(mol);
  get_element_type(mol,1,the_type);
  fprintf(file1," # TITLE\n");
  fprintf(file1,"%3d  0DIST  0  0  0\n",Atoms);
  fprintf(file1,"  EL\n");
  fprintf(file1,"0.,0.,0., %s\n",the_type);
  for (i = 2; i <= Atoms; i++)
  {
    get_element_type(mol,i,the_type);
    if (T(i) < 0.0)
      T(i) += 360.0;
    fprintf(file1,"%2d,%d,%2s%7.3f,%7.3f,%7.3f\n",
           NA(i),i,the_type,R(i),W(i),T(i));
  }
  return(TRUE);
}

void add_dummy_atoms(ums_type *mol)
{
  int i;

  Atoms += 2;
  reinitialize_ums(&mol);

  i = Atoms-1;
  X(i) = 0.0;
  Y(i) = 0.0;
  Z(i) = 1.0;
  i = Atoms;
  X(i) = 1.0;
  Y(i) = 0.0;
  Z(i) = 0.0;
  Atoms -= 2;
}

void center_at_first_atom(ums_type *mol)
{
  int i;
  double tempX, tempY, tempZ;

  tempX = X(1);
  tempY = Y(1);
  tempZ = Z(1);

  for (i = 1; i <= Atoms; i++)
  {
    X(i) = X(i) - tempX;
    Y(i) = Y(i) - tempY;
    Z(i) = Z(i) - tempZ;
  }
}

void set_hilderbrandt_connections(ums_type *mol)
{
  int i,j,k = 0;
  double sum,r;

  NA(1) = Atoms+1;
  NB(1) = Atoms+2;
  NC(1) = 0;
  NB(2) = Atoms+1;
  NC(2) = Atoms+2;
  NC(3) = Atoms+1;

  for (i = 2; i <= Atoms; i++)
  {
    sum = 100.0;
    for (j = 1; j < i; j++)
    {
      r = distance(Point(i),Point(j));
      if ((r < sum) && (NA(j) != j) && (NB(j) != j))
      {
        sum = r;
        k = j;
      }
    }
    NA(i) = k;
  }
  for (i = 3; i <= Atoms; i++)
  {
    NB(i) = NA(NA(i));
  }
  for (i = 4; i <= Atoms; i++)
  {
    NC(i) = NB(NB(i));
  }
}

void set_hilderbrandt_geometry(ums_type *mol)
{
  int i;

  R(1) = 0.0;
  W(1) = 0.0;
  T(1) = 0.0;
  W(2) = 0.0;
  T(2) = 0.0;
  T(3) = 0.0;

  for (i = 2; i <= Atoms; i++)
  {
    W(i) = bond_angle(Point(i),Point(NA(i)),Point(NB(i)));
    T(i) = torsion(Point(i),Point(NA(i)),Point(NB(i)),Point(NC(i)));
    R(i) = distance(Point(i),Point(NA(i)));
  }
}
