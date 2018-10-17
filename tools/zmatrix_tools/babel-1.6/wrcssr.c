/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
----------------------------------------------------------------------
FILE : wrcssr.c
AUTHOR(S) : Pat Walters
DATE : 2-94
PURPOSE : Routines to write a CSD CSSR type file
******/

#include "bbltyp.h"

int 
write_cssr(FILE *file1, ums_type *mol)
{ 
  int i,j;
  char type_name[5];
  pdb_type_rec *pdb_types;

  int result;

  pdb_types = (pdb_type_rec *)malloc((Atoms + 1) * sizeof(pdb_type_rec));
  
  /*
  coord_type dim;
  dim = calc_cell_dimensions(mol); 
  make_fractional(mol,dim);
  */
  if (mol->fract == NULL)
  {
    fprintf(file1,
	    " REFERENCE STRUCTURE = 00000   A,B,C =  %6.3f  %6.3f  %6.3f\n",
	    1.0,1.0,1.0);
    fprintf(file1,
	    "   ALPHA,BETA,GAMMA =  90.000  90.000  90.000    SPGR =    P1\n");
    fprintf(file1,"%4d\n\n",Atoms);
  }
  else
  {
    fprintf(file1,
	    " REFERENCE STRUCTURE = 00000   A,B,C =  %6.3f  %6.3f  %6.3f\n",
	    mol->fract->A,mol->fract->B,mol->fract->C);
    fprintf(file1,
	    "   ALPHA,BETA,GAMMA =  %6.4f  %6.3f  %6.3f    SPGR =    ??\n",
	    mol->fract->Alpha,mol->fract->Beta,mol->fract->Gamma);
    fprintf(file1,"%4d\n\n",Atoms);
  }
  
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    strcpy(pdb_types[i].name,type_name);
    assign_pdb_number(pdb_types,i);

    fprintf(file1," %3d%2s%-3d  %8.5f  %8.5f  %8.5f ",
	    i,
	    type_name,
	    pdb_types[i].number,
	    X(i),
	    Y(i),
	    Z(i));
    for (j = 0; j < Valence(i); j++)
      fprintf(file1,"%4d",Connection(i,j));
    fprintf(file1,"\n");
  }
  return(TRUE);
}


coord_type calc_cell_dimensions(ums_type *mol)
{
  int i;
  double Xmax = -999999.0;
  double Ymax = -999999.0;
  double Zmax = -999999.0;
  double Xmin = 999999.0;
  double Ymin = 999999.0;
  double Zmin = 999999.0;
  coord_type point;

  for (i = 1; i < Atoms; i++)
  {
    if (X(i) > Xmax) Xmax = X(i);
    if (Y(i) > Ymax) Ymax = Y(i);
    if (Z(i) > Zmax) Zmax = Z(i);
    if (X(i) < Xmin) Xmin = X(i);
    if (Y(i) < Ymin) Ymin = Y(i);
    if (Z(i) < Zmin) Zmin = Z(i);
  }
  point.x = Xmax - Xmin;
  point.y = Ymax - Ymin;
  point.z = Zmax - Zmin;

  return(point);
}

    
void make_fractional(ums_type *mol, coord_type dim)
{
  int i;
  
  for (i = 1; i <= Atoms; i++)
  {
    X(i) = X(i)/dim.x;
    Y(i) = Y(i)/dim.y;
    Z(i) = Z(i)/dim.z;
  }
}









