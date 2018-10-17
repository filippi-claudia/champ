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

FILE : rdfract.c
AUTHOR(S) : Pat Walters
DATE : 6-24-93
PURPOSE : Routines to read a CSD CSSR file

******/

#include "bbltyp.h"

#define DELIMS " \t\n"

int 
  read_csd_fractional(FILE *file1, ums_type *mol)
{
  char xyz_line[BUFF_SIZE];
  int i,j;
  int result;
  int tokens;
  matrix_3x3 m;
  fract_type f;
  
  fgets(xyz_line,sizeof(xyz_line), file1);
  sscanf(&xyz_line[39],"%lf%lf%lf",&f.A,&f.B,&f.C);
  fgets(xyz_line,sizeof(xyz_line), file1);
  sscanf(&xyz_line[21],"%lf%lf%lf",&f.Alpha,&f.Beta,&f.Gamma);
  fgets(xyz_line,sizeof(xyz_line), file1);
  sscanf(xyz_line,"%d",&Atoms);
  fill_orth_matrix(&f,&m);
  ShowProgress(Atoms,"Reading Atoms");
  result = initialize_ums(&mol);
  fgets(xyz_line,sizeof(xyz_line), file1);
  strncpy(Title,xyz_line,sizeof(Title));
  strip_return(Title);
  for (i = MIN_ATOM; i <= Atoms;i ++)
  {
    UpdateProgress();
    fgets(xyz_line,sizeof(xyz_line), file1);
    sscanf(xyz_line,"%*s %s %lf %lf %lf",
	   Type(i),
	   &X(i),
	   &Y(i),
	   &Z(i));
    clean_atom_type(Type(i));
    fract_to_cart(&Point(i),&m);
    tokens = count_tokens(xyz_line,DELIMS);
    if (tokens >= 13)
    {
      for (j = 0; j < 8; j++)
      {
	sscanf(&xyz_line[41 + j * 4],"%d",&Connection(i,j));
	if (Connection(i,j) != 0)
	  Valence(i)++;
      }
    }
    else
    {
      Valence(i) = tokens - 5;
      for (j = 0; j < (Valence(i)); j++)
      {
	sscanf(&xyz_line[41 + j * 4],"%d",&Connection(i,j));
	if (Connection(i,j) == 0)
	{
	  Valence(i) = 0;
	  break;
	}
      }
    }
  }
  
  result = assign_types(mol);
  result = build_connection_table(mol);
  return(TRUE);
}


int 
  read_fform_fract(FILE *file1, ums_type *mol)
{
  char xyz_line[BUFF_SIZE];
  int i;
  int result;
  fract_type f;
  matrix_3x3 m;
  
  fgets(xyz_line,sizeof(xyz_line), file1);
  sscanf(xyz_line,"%d%lf%lf%lf%lf%lf%lf",&Atoms,&f.Alpha,&f.Beta,&f.Gamma,&f.A,&f.B,&f.C);
  fill_orth_matrix(&f,&m);
  result = initialize_ums(&mol);
  fgets(xyz_line,sizeof(xyz_line), file1);
  for (i = MIN_ATOM; i <= Atoms;i ++)
  {
    fgets(xyz_line,sizeof(xyz_line), file1);
    sscanf(xyz_line,"%s %lf %lf %lf",
	   Type(i),
	   &X(i),
	   &Y(i),
	   &Z(i));
    clean_atom_type(Type(i));
    fract_to_cart(&Point(i),&m);
  }
  
  result = assign_radii(mol);
  result = assign_bonds(mol);
  result = assign_types(mol);
  result = build_connection_table(mol);
  assign_bond_order(mol);
  
  return(TRUE);
}




/***
  Fill in the elements of the orthogonalization matrix
  AUTHOR : Pat Walters
  DATE : 6-24-93
  void fill_orth_matrix(fract_type *f, matrix_3x3 *m)
  f - pointer to structure with A,B,C,Alpha,Beta,Gamma
  m - pointer to the orthogonaliztion matrix
  
  ***/

void fill_orth_matrix(fract_type *f, matrix_3x3 *m)
{
  double V;
  double Alpha,Beta,Gamma;
  
  Alpha = f->Alpha * DEG_TO_RAD;
  Beta = f->Beta * DEG_TO_RAD;
  Gamma = f->Gamma * DEG_TO_RAD;
  
  V= 1 - SQUARE(cos(Alpha)) - SQUARE(cos(Beta)) - SQUARE(cos(Gamma)) 
    + 2 * cos(Alpha) * cos(Beta) *  cos(Gamma);
  V = sqrt(V)/sin(Gamma);
  
  m->a1 = f->A;
  m->b1 = f->B * cos(Gamma);
  m->c1 = f->C * cos(Beta);
  m->a2 = 0.0;
  m->b2 = f->B * sin(Gamma);
  m->c2 = f->C * (cos(Alpha)-cos(Beta)*cos(Gamma))/sin(Gamma);
  m->a3 = 0.0;
  m->b3 = 0.0;
  m->c3 = f->C * V;
}


/***
  Convert a point in fractional coordinates to a point in
  cartesian coordinates
  AUTHOR : Pat Walters
  DATE : 6-24-93
  
  void fract_to_cart(coord_type *p, matrix_3x3 *m)
  p - pointer to the coordinate structure
  m - the ortogonalization matrix
  ***/


void fract_to_cart(coord_type *p, matrix_3x3 *m)
{
  vect_type v;
  
  v.x = p->x;
  v.y = p->y;
  v.z = p->z;
  
  mat_3x3_dot_vect(m, &v);
  
  p->x = v.x;
  p->y = v.y;
  p->z = v.z;
}










