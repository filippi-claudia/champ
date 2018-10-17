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

FILE : assbnd.c
AUTHOR(S) : Pat Walters
DATE : 10-15-92
PURPOSE : Routines to determine which atoms are bonded based on 
a comparison of bond distances and covalent radii
HISTORY : 
8-5-93 Modified assign_bonds so that it uses as coordinate presort
as suggested by Anders Sundin  
******/
#include "bbltyp.h"

/************************************************************************
look up the radius and maximum # of connections in cov_radius and max_bonds
*************************************************************************/
extern element_type *elements; 
static warning wstr;

/*#define SLOP_FACTOR 0.2*/
#define SLOP_FACTOR 0.45
#define CLOSE_CUTOFF 2.0

#undef DEBUG

int 
assign_radii(ums_type *mol)
{
  int i, j;
  int found;

  for (i = MIN_ATOM; i <= Atoms; i++)
  {
    found = FALSE;
    for (j = 0; j < MAX_ELEMENTS; j++)
    {
      if (strncmp(Type(i),elements[j].name,2) == 0) 
      {
	Radius(i) = elements[j].cov_rad;
	BORadius(i) = elements[j].bond_ord_rad;
	Max_bonds(i) = elements[j].max_bonds;
	found = TRUE;
	break;
      }
    }
    if (found == FALSE) 
    {
      sprintf(wstr,"No bond length information for atom %d, type = %s",
	     i,Type(i));
      show_warning(wstr);
    }
  }
  return(1);
}

void fast_radii(ums_type *mol)
{
  int i;
  
  for (i = 1; i <= Atoms; i++)
  {
    Radius(i) = elements[Atomic_number(i)].cov_rad;
  }
}


int 
assign_hybrid_radii(ums_type *mol)
{
  int i,j;
  int found;
  char atm_type[5];

  for (i = MIN_ATOM; i <= Atoms; i++)
  {
    get_element_type(mol,i,atm_type);
    clean_atom_type(atm_type);
    found = FALSE;
    for (j = 0; j < MAX_ELEMENTS; j++)
    {
      if (EQn(atm_type,elements[j].name,2))
      {
	Radius(i) = elements[j].cov_rad;
	BORadius(i) = elements[j].bond_ord_rad;
	Max_bonds(i) = elements[j].max_bonds;
	found = TRUE;
	break;
      }
    }
    if (!found) 
    {
      sprintf(wstr,"No radius information for atom %d, type = %s",
	     i,atm_type);
      show_warning(wstr);
    }
  }
  return(1);
}


int assign_bonds(ums_type *mol)
{
  int i,j;
  temp_atom_rec *temp;
  double xDiff, yDiff, zDiff, MaxLength, dist2;

  Bonds = 0;

  temp = (temp_atom_rec *)malloc(Atoms * sizeof(temp_atom_rec));

  if (!temp)
    fatal_error("Memory allocation error");

  for (i = 0; i < Atoms; i++)
  {
    j = i + 1;
    temp[i].x = X(j);
    temp[i].y = Y(j);
    temp[i].z = Z(j);
    temp[i].num = j;
  }
  
  qsort(temp,Atoms,sizeof(temp_atom_rec),QSORT_PROTO zsort_atoms);
  
  for (i = 0; i < Atoms; i++)
    for (j = 0; j < Atoms; j++)
    {
      if (i != j)
      {
	zDiff = temp[j].z - temp[i].z;
	MaxLength = Radius(temp[i].num) + Radius(temp[j].num) + SLOP_FACTOR;
/*	if (zDiff > MaxLength)
	  break;
*/
	xDiff = temp[i].x - temp[j].x;
	yDiff = temp[i].y - temp[j].y;
	dist2 = SQUARE(xDiff) + SQUARE(yDiff) + SQUARE(zDiff);
	if (dist2 < SQUARE(MaxLength)) 
	{
	  if (temp[i].num < temp[j].num)
	  {
	    Start(Bonds) = temp[i].num;
	    End(Bonds) = temp[j].num;
	    Bonds ++;
	  }
	}
      }
    }
  qsort(mol->connections,Bonds,sizeof(connect_type),QSORT_PROTO sort_connections);
  for (i = 1; i <= Atoms; i++)
    Valence(i) = 0;
  dissect_connection_table(mol);
  check_bonds(mol);
  free(temp); 
  return(TRUE);
}

int assign_pdb_bonds(ums_type *mol)
{
  int i,j;
  double xDiff, yDiff, zDiff, MaxLength, dist2;
  int diff;
  int res_dist;

  for (i = 1; i < Atoms; i++)
  {
    for (j = i+1; j <= Atoms; j++)
    {
      if (ChainNum(i) == ChainNum(j))
      {
        if (!bonded(mol, i, j))
        {
          if (Type(i)[0] == 'H')
  	    res_dist = 1;
          else
	    res_dist = 2;
	  
	  if (ResNum(j) > ResNum(i) + res_dist) break;
	  
          diff = abs(ResNum(i) - ResNum(j));
          if (diff < res_dist)
          {
	    MaxLength = Radius(i) + Radius(j) + SLOP_FACTOR;
	    xDiff = X(i) - X(j);
	    yDiff = Y(i) - Y(j);
	    zDiff = Z(i) - Z(j);
	    dist2 = SQUARE(xDiff) + SQUARE(yDiff) + SQUARE(zDiff);
	    if (dist2 < SQUARE(MaxLength)) 
	    {
              Connection(i,Valence(i)) = j;
              Connection(j,Valence(j)) = i;
	      Valence(i)++;
	      Valence(j)++;
	      Start(Bonds) = i;
	      End(Bonds) = j;
	      Bonds ++;
	    }
          }
        }
      }
    }
  }

  return(TRUE);
}


void estimate_bond_order(ums_type *mol)
{
  double ratio,dist,cov_sum;
  int bo;
  int i;
  
  for (i = 0; i < Bonds; i++)
  {
    bo = 1;
    dist = distance(Point(Start(i)),Point(End(i)));
    cov_sum = BORadius(Start(i)) + BORadius(End(i));
    ratio = dist/cov_sum;
    if (ratio <= 0.81)
      bo = 3;
    else
      if (ratio <= 0.94)
	bo = 2;
#ifdef DEBUG
    printf("i = %d Start = %d End = %d r1 = %10.3f r2 = %10.3f dist = %10.2f sum = %10.2f ratio = %10.2f bond order = %4d\n",
	   i,Start(i),End(i),BORadius(Start(i)),BORadius(End(i)),dist,cov_sum,ratio,bo);
#endif
    Bond_order(i) = bo;
  }
}


void check_bonds(ums_type *mol)
{
  int i;
  
  for (i = 1; i <= Atoms; i++)
  {
    if (Valence(i) > Max_bonds(i))
    {
      sprintf(wstr,
	      "Valence > maximum for atom=%d type=%s valence=%d max=%d",
	      i,
	      Type(i),
	      Valence(i),
	      Max_bonds(i));
      print_bad_connections(mol,i);
      show_warning(wstr);
    }
    if (Valence(i) == 0)
    {
	sprintf(wstr,"atom %d is isolated, type=%s",i,Type(i));
	if (find_close_atoms(mol,i) == FALSE) 
	  show_warning(wstr);
    }
  }
}


int find_close_atoms(ums_type *mol, int num)
{
  int i;
  double dist;
  int found = FALSE;
  
  for (i = 1; i <= Atoms; i++)
  {
    dist = distance(Point(i),Point(num));
    if (i != num)
    {
      if (dist < (Radius(i) + Radius(num) + SLOP_FACTOR)) 
      {
	found = TRUE;
	Connection(num,Valence(num)) = i;
	Valence(num)++;
	Connection(i,Valence(i)) = num;
	Valence(i)++;
      }
    }
  }
  return(found);
}



int zsort_atoms(temp_atom_rec *a, temp_atom_rec *b)
{
  if (a->z > b->z) 
    return(1);
  else 
    if (a->z < b->z)
      return(-1);
  else
    return(0);
}

int sort_connections(connect_type *a, connect_type *b)
{
  if (a->start > b->start)
    return(1);
  else 
    if (a->start < b->start)
      return(-1);
  else
    if (a -> end > b->end)
      return(1);
  else
    return(-1);
}

double 
distance(coord_type first,coord_type second)
{
  double dist;
  dist = sqrt(SQUARE(first.x - second.x) + SQUARE(first.y - second.y) + SQUARE(first.z - second.z));
  return(dist);
}


int is_element(char *ele)
{
  int i;
  clean_atom_type(ele);
  for (i = 0; i < MAX_ELEMENTS; i++)  
  {
    if (EQ(ele,elements[i].name))
      return(TRUE);
  }
  return(FALSE);
}

int 
get_atomic_number(char *sym)
{
  int i;

  clean_atom_type(sym);
  for (i = 0; i < MAX_ELEMENTS; i++)
  {
    if (EQ(elements[i].name,sym))
      return(i);
  }
  sprintf(wstr,"Could not find atomic number for %s",sym);
  show_warning(wstr);
  return(0);
}













