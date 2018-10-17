/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
----------------------------------------------------------------------------

FILE : orient.c
AUTHOR(S) : Pat Walters
DATE : 9-1-94
PURPOSE : routine to orient a structure along the x and y axes
******/


#include "bbltyp.h"
#define DELIMS " ,\t\n"
#define Primary    1
#define Secondary  2
#define Tertiary   3

void set_std_orientation(ums_type *mol)
{
  center_at_atom(mol,2);
  orient_ums(mol,1,2,1,3);
}


void get_orientation_matrix(ums_type *mol, int x1, int x2, int y1, int y2, matrix_3x3 *m)
{
  vect_type r1,r2,v1,v2,v3;
  
  pts_2_vect(mol,&r1,x1,x2);
  pts_2_vect(mol,&r2,y1,y2);
  cross_prod(&r1,&r2,&v3);
  normalize_vect(&v3);
  cross_prod(&v3,&r1,&v2);
  normalize_vect(&v2);
  cross_prod(&v2,&v3,&v1);
  normalize_vect(&v1);
  
  m->a1 = v1.x; m->b1 = v1.y; m->c1 = v1.z;
  m->a2 = v2.x; m->b2 = v2.y; m->c2 = v2.z;
  m->a3 = v3.x; m->b3 = v3.y; m->c3 = v3.z;  
}

void ums_plus_vector(ums_type *mol, vect_type *v)
{
  int i;
  
  for (i = 1; i <= Atoms; i++)
  {
    X(i) += v->x;
    Y(i) += v->y;
    Z(i) += v->z;
  }
}



void ums_dot_matrix(ums_type *mol, matrix_3x3 *m)
{
  int i;
  vect_type temp;
  
  for (i = 1; i <= Atoms; i++)
  {
    coord_to_vector(&temp,&Point(i));
    mat_3x3_dot_vect(m,&temp);
    vector_to_coord(&Point(i),&temp);
  }
}


void orient_ums(ums_type *mol, int x1, int x2, int y1, int y2)
{
  int i;
  vect_type r1,r2,v1,v2,v3, temp;
  matrix_3x3 m;
  
  pts_2_vect(mol,&r1,x1,x2);
  pts_2_vect(mol,&r2,y1,y2);
  cross_prod(&r1,&r2,&v3);
  normalize_vect(&v3);
  cross_prod(&v3,&r1,&v2);
  normalize_vect(&v2);
  cross_prod(&v2,&v3,&v1);
  normalize_vect(&v1);

  m.a1 = v1.x; m.b1 = v1.y; m.c1 = v1.z;
  m.a2 = v2.x; m.b2 = v2.y; m.c2 = v2.z;
  m.a3 = v3.x; m.b3 = v3.y; m.c3 = v3.z;

  for (i = 1; i <= Atoms; i++)
  {
    coord_to_vector(&temp,&Point(i));
    mat_3x3_dot_vect(&m,&temp);
    vector_to_coord(&Point(i),&temp);
  }
}  

void coord_to_vector(vect_type *v, coord_type *c)
{
  v->x = c->x;
  v->y = c->y;
  v->z = c->z;
}

void vector_to_coord(coord_type *c, vect_type *v)
{
  c->x = v->x;
  c->y = v->y;
  c->z = v->z;
}

void center_at_atom(ums_type *mol, int atom)
{
  int i;
  double tempX, tempY, tempZ;
  
  tempX = X(atom);
  tempY = Y(atom);
  tempZ = Z(atom);
  
  for (i = 1; i <= Atoms; i++)
  {
    X(i) -= tempX;
    Y(i) -= tempY;
    Z(i) -= tempZ;
  }
}

void center_at_atoms(ums_type *mol,int atom[],int count)
{
  int i;
  double tempX, tempY, tempZ;

  tempX = tempY = tempZ = 0.0;
  for (i = 0;i < count;i++)
  {
    tempX += X(atom[i]);
    tempY += Y(atom[i]);
    tempZ += Z(atom[i]);
  }

  tempX /= (double)count;
  tempY /= (double)count;
  tempZ /= (double)count;
  
  for (i = 1; i <= Atoms; i++)
  {
    X(i) -= tempX;
    Y(i) -= tempY;
    Z(i) -= tempZ;
  }
}

void center_at_origin(ums_type *mol, vect_type *v)
{
  int i;
  double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
  
  for (i = 1; i <= Atoms; i++)
  {
    sumX += X(i);
    sumY += Y(i);
    sumZ += Z(i);
  }
  
  sumX /= (double)Atoms;
  sumY /= (double)Atoms;
  sumZ /= (double)Atoms;
  
  for (i = 1; i <= Atoms; i++)
  {
    X(i) -= sumX;
    Y(i) -= sumY;
    Z(i) -= sumZ;
  }

  v->x = sumX;
  v->y = sumY;
  v->z = sumZ;
}

void AlignMol(ums_type *mol)
{
  int i;
  int *cList,*pList,*sList;
  int cCount,primCount,secCount;
  vect_type r1,r2,v1,v2,v3;
  int x,y,z;
  double primary,secondary,tertiary;
  char buffer[100];
  char axis;

  fprintf(stdout,"Input the atom(s) to be placed at the origin\n");
  fgets(buffer,sizeof(buffer),stdin);
  cCount = count_tokens(buffer,DELIMS);
  if (cCount)
  {
    cList = (int *)malloc(sizeof(int) * cCount);
    if (!cList)
      fatal_error("Unable to allocate memory for array");
    
    for (i = 0;i < cCount;i++)
      cList[i] = atoi(gettoken(buffer,DELIMS,i+1));
  }
  else
    fatal_error("Atoms to be placed at the origin must be defined");

  x = y = z = 0;

  fprintf(stdout,"What is the primary alignment axis (x, y, or z)?\n");
  fgets(buffer,sizeof(buffer),stdin);
  lowercase(buffer);
  axis = '\0';
  sscanf(buffer,"%c",&axis);

  if (!axis)
    fatal_error("You must input x, y, or z!");

  switch (axis)
  {
  case 'x':
    x = Primary;
    break;
  case 'y':
    y = Primary;
    break;    
  case 'z':
    z = Primary;
    break;
 }

  fprintf(stdout,"Enter the atoms to be aligned along the primary axis\n");
  fgets(buffer,sizeof(buffer),stdin);
  primCount = count_tokens(buffer,DELIMS);
  pList = (int *)malloc(sizeof(int) * primCount);
  if (!pList)
    fatal_error("Unable to allocate memory for array");
  
  for (i = 0;i < primCount;i++)
    pList[i] = atoi(gettoken(buffer,DELIMS,i+1));
  

  fprintf(stdout,"What is the secondary alignment axis (x, y, or z)?\n");
  fgets(buffer,sizeof(buffer),stdin);
  lowercase(buffer);
  axis = '\0';
  sscanf(buffer,"%c",&axis);

  if (!axis)
    fatal_error("You must input x, y, or z!");

  switch (axis)
  {
  case 'x':
    if (x)
      fatal_error("X has already been defined");
    x = Secondary;
    break;
  case 'y':
    if (y)
      fatal_error("Y has already been defined");
    y = Secondary;
    break;    
  case 'z':
    if (z)
      fatal_error("Z has already been defined");
    z = Secondary;
    break;
 }

  fprintf(stdout,"Enter the atoms to be aligned along the secondary axis\n");
  fgets(buffer,sizeof(buffer),stdin);
  secCount = count_tokens(buffer,DELIMS);
  sList = (int *)malloc(sizeof(int) * primCount);
  if (!sList)
    fatal_error("Unable to allocate memory for array");
  
  for (i = 0;i < secCount;i++)
    sList[i] = atoi(gettoken(buffer,DELIMS,i+1));

  
  if (!(x == Primary || y == Primary || z == Primary))
    fatal_error("A primary axis must be defined");

  if (!(x == Secondary || y == Secondary || z == Secondary))
    fatal_error("A secondary axis must be defined");
  
  if (!x)
    x = Tertiary;
  if (!y)
    y = Tertiary;
  if (!z)
    z = Tertiary;

/*  pts_2_vect(mol,&r1,x1,x2);
  pts_2_vect(mol,&r2,y1,y2);*/  

  center_at_atoms(mol,cList,cCount);

  r1.x = r1.y = r1.z = 0;
  for (i = 0;i < primCount;i++)  
  {
    r1.x += X(pList[i]);
    r1.y += Y(pList[i]);
    r1.z += Z(pList[i]);
  }
  r1.x /= (double) primCount;
  r1.y /= (double) primCount;
  r1.z /= (double) primCount;

  r2.x = r2.y = r2.z = 0;
  for (i = 0;i < secCount;i++)  
  {
    r2.x += X(sList[i]);
    r2.y += Y(sList[i]);
    r2.z += Z(sList[i]);
  }
  r2.x /= (double) secCount;
  r2.y /= (double) secCount;
  r2.z /= (double) secCount;


  cross_prod(&r1,&r2,&v3);
  normalize_vect(&v3);
  cross_prod(&v3,&r1,&v2);
  normalize_vect(&v2);
  cross_prod(&v2,&v3,&v1);
  normalize_vect(&v1);

  for (i = 1; i <= Atoms; i++)
  {
    primary = v1.x * X(i) + v1.y * Y(i) + v1.z * Z(i);
    secondary = v2.x * X(i) + v2.y * Y(i) + v2.z * Z(i);
    tertiary = v3.x * X(i) + v3.y * Y(i) + v3.z * Z(i);

    switch (x)
    {
    case Primary:
      X(i) = primary;
      break;
    case Secondary:
      X(i) = secondary;
      break;
    case Tertiary:
      X(i) = tertiary;
      break;
    }

    switch (y)
    {
    case Primary:
      Y(i) = primary;
      break;
    case Secondary:
      Y(i) = secondary;
      break;
    case Tertiary:
      Y(i) = tertiary;
      break;
    }

    switch (z)
    {
    case Primary:
      Z(i) = primary;
      break;
    case Secondary:
      Z(i) = secondary;
      break;
    case Tertiary:
      Z(i) = tertiary;
      break;
    }
  }
}

void find_mol_center(ums_type *mol, vect_type *center, int use_wts)
{
  int i;
  double wnorm = 0.0;
  double maxX = -99e99, maxY = -99e99, maxZ = -99e99;
  double minX =  99e99, minY =  99e99, minZ =  99e99;

  center->x = 0.0;
  center->y = 0.0;
  center->z = 0.0;

  if (use_wts)
    {
      get_atomic_weights(mol);
    }
  else
    {
      for (i = 1; i <= Atoms; i++)
	{
	  Double(i) = 1.0;
	}
    }
  
  for (i = 1; i <= Atoms; i++) 
    {
      center->x += X(i) * Double(i);
      center->y += Y(i) * Double(i);
      center->z += Z(i) * Double(i);
      wnorm += Double(i);
    }

  center->x /= wnorm;
  center->y /= wnorm;
  center->z /= wnorm;
}









