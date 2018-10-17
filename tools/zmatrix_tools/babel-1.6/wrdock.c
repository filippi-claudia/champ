/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
------------------------------------------------------------------------------

FILE : wrxyz.c
AUTHOR(S) : Pat Walters
DATE : 11-92
PURPOSE : Routines to write a dock database file
******/

#include "bbltyp.h"

static int first_one = TRUE;

#define DOCK_CHG_FACTOR 1000.0

int 
write_dock(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  int result;
  int num, heavy_count, h_count;
  double Xmin = 999999.0;
  double Ymin = 999999.0;
  double Zmin = 999999.0;
  vect_type v;
  ums_type *temp = NULL;
  char temp_title[BUFF_SIZE];

  if (first_one)
  {
    fprintf(file1,"%s\n","DOCK 3.5 ligand_atoms");
    first_one = FALSE;
  }

  type_hydrogens(mol);

  for(i = 1;i <= Atoms; i++)
  {
    if (X(i) < Xmin) Xmin = X(i);
    if (Y(i) < Ymin) Ymin = Y(i);
    if (Z(i) < Zmin) Zmin = Z(i);
  }

  v.x = -Xmin;
  v.y = -Ymin;
  v.z = -Zmin;
  
  ums_plus_vector(mol,&v);

  heavy_count = total_heavy_atoms(mol);
  h_count = Atoms - heavy_count;
  
  fprintf(file1,"%-51s%-51s\n",Title,Title);
  fprintf(file1,"%3d%3d%3d%12.4f%6d%9.2f\n",
	 Atoms,heavy_count,h_count,0.0,0,0.0);

  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"DOK",Type(i),type_name,zero);
    fprintf(file1,"%5.0f%5.0f%5.0f%2d%5.0f%1d%3d\n",
	    X(i)*1000.0,
	    Y(i)*1000.0,
	    Z(i)*1000.0,
	    atoi(type_name),
	    Charge(i)*DOCK_CHG_FACTOR,
	    0,0);
  }
  return(TRUE);
}


int total_heavy_atoms(ums_type *mol)
{
  int i;
  int heavy = 0;
  
  for (i = 1; i <= Atoms; i++)
  {
    if (Atomic_number(i) != 1)
      heavy++;
  }

  return(heavy);
}


      













