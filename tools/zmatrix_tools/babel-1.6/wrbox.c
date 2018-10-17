/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
------------------------------------------------------------------------------

FILE : wrbox.c
AUTHOR(S) : Pat Walters
DATE : 10-31-96  Boo !!
PURPOSE : Routines to write a pdb box for dock 3.5
******/

#include "bbltyp.h"

int 
write_box(FILE *fp, ums_type *mol)
{ 
  vect_type center,min,max,mid,dim;
  double margin;
  
  /*  
  if (HasAttribute(ControlMem) && (InputKeywords))
    margin = atof(InputKeywords);
  else
  */
  margin = 10.0;
  
  find_mol_center(mol,&center,FALSE);
  find_extents(mol,&min,&max);

  min.x -= margin;
  min.y -= margin;
  min.z -= margin;
  max.x += margin;
  max.y += margin;
  max.z += margin;

  dim.x = max.x - min.x;
  dim.y = max.y - min.y;
  dim.z = max.z - min.z;

  mid.x = (max.x + min.x)/2.0;
  mid.y = (max.y + min.y)/2.0;
  mid.z = (max.z + min.z)/2.0;

  fprintf(fp,"HEADER    CORNERS OF BOX\n");
  fprintf(fp,"REMARK    CENTER (X Y Z)      %10.3f %10.3f %10.3f\n",
	  mid.x,mid.y,mid.z);
  fprintf(fp,"REMARK    DIMENSIONS (X Y Z)  %10.3f %10.3f %10.3f\n",
	  dim.x,dim.y,dim.z);
  fprintf(fp,"ATOM      1  DUA BOX     1    %8.3f%8.3f%8.3f\n",
	  mid.x - dim.x/2.0, mid.y - dim.y/2.0, mid.z - dim.z/2.0);
  fprintf(fp,"ATOM      2  DUB BOX     1    %8.3f%8.3f%8.3f\n",
	  mid.x + dim.x/2.0, mid.y - dim.y/2.0, mid.z - dim.z/2.0);
  fprintf(fp,"ATOM      3  DUC BOX     1    %8.3f%8.3f%8.3f\n",
	  mid.x + dim.x/2.0, mid.y - dim.y/2.0, mid.z + dim.z/2.0);
  fprintf(fp,"ATOM      4  DUC BOX     1    %8.3f%8.3f%8.3f\n",
	  mid.x - dim.x/2.0, mid.y - dim.y/2.0, mid.z + dim.z/2.0);
  fprintf(fp,"ATOM      5  DUD BOX     1    %8.3f%8.3f%8.3f\n",
	  mid.x - dim.x/2.0, mid.y + dim.y/2.0, mid.z - dim.z/2.0);
  fprintf(fp,"ATOM      6  DUE BOX     1    %8.3f%8.3f%8.3f\n",
	  mid.x + dim.x/2.0, mid.y + dim.y/2.0, mid.z - dim.z/2.0);
  fprintf(fp,"ATOM      7  DUF BOX     1    %8.3f%8.3f%8.3f\n",
	  mid.x + dim.x/2.0, mid.y + dim.y/2.0, mid.z + dim.z/2.0);
  fprintf(fp,"ATOM      8  DUG BOX     1    %8.3f%8.3f%8.3f\n",
	  mid.x - dim.x/2.0, mid.y + dim.y/2.0, mid.z + dim.z/2.0);

  fprintf(fp,"CONECT    1    2    4    5\n");
  fprintf(fp,"CONECT    2    1    3    6\n");
  fprintf(fp,"CONECT    3    2    4    7\n");
  fprintf(fp,"CONECT    4    1    3    8\n");
  fprintf(fp,"CONECT    5    1    6    8\n");
  fprintf(fp,"CONECT    6    2    5    7\n");
  fprintf(fp,"CONECT    7    3    6    8\n");
  fprintf(fp,"CONECT    8    4    5    7\n");
  
  return(TRUE);
}


void find_extents(ums_type *mol, vect_type *min, vect_type *max)
{
  int i;
  
  max->x = max->y = max->z = -99999.0;
  min->x = min->y = min->z =  99999.0;

  for (i = 1; i <= Atoms; i++)
    {
      min->x = MIN(X(i),min->x);
      min->y = MIN(Y(i),min->y);
      min->z = MIN(Z(i),min->z);
      max->x = MAX(X(i),max->x);
      max->y = MAX(Y(i),max->y);
      max->z = MAX(Z(i),max->z);
    }
}

  
      



