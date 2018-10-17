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

FILE : rdmolen.c
AUTHOR(S) : Pat Walters
DATE : 4-93
PURPOSE : routines to read file produced by the MOLIN program
******/


#include "bbltyp.h"

int 
read_molin(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i=0;
  int result;
  matrix_3x3 m;
  int keep_fractional = FALSE;
  coord_type *tmp;
  
  uppercase(InputKeywords);
  if EQ(InputKeywords,"NOCART")
    keep_fractional = TRUE;
  while(fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (count_tokens(the_line,"\n\t ") == 4)
    {
      i++;
    }
  }
  Atoms = i ;
  tmp = (coord_type *)malloc((Atoms + 1) * sizeof(coord_type));
  result = initialize_ums(&mol);
  initialize_fractional(&mol);
  rewind(file1);
  fgets(the_line,sizeof(the_line), file1);
  if ((count_tokens(the_line,"\t\n ") == 6))
    sscanf(the_line,"%lf%lf%lf%lf%lf%lf",
	   &mol->fract->A,&mol->fract->B,&mol->fract->C,
	   &mol->fract->Alpha,&mol->fract->Beta,&mol->fract->Gamma);
  else
  {
    get_cell_params(mol->fract);
    rewind(file1);
  }
  fill_orth_matrix(mol->fract,&m);
  ShowProgress(Atoms,"Reading Atoms");
  i = 1;
  while(fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (count_tokens(the_line,"\t\n ") == 4)
    {
      UpdateProgress();
      sscanf(the_line,"%s %lf %lf %lf",
	     Type(i),
	     &X(i),
	     &Y(i),
	     &Z(i));      
      clean_atom_type(Type(i)); 
      tmp[i].x = X(i);
      tmp[i].y = Y(i);
      tmp[i].z = Z(i);
      fract_to_cart(&Point(i),&m);
      i++;
      }
    }
  result = assign_radii(mol);
  result = assign_bonds(mol);
  result = assign_types(mol);
  assign_bond_order(mol);
  if (keep_fractional)
  {
    for (i = 1; i <= Atoms; i++)
    {
      X(i) = tmp[i].x;
      Y(i) = tmp[i].y;
      Z(i) = tmp[i].z;
    }
  }
  else
    if (mol->fract)
    	free(mol->fract);
  if (tmp)
    free(tmp);
  return(TRUE);
}


void get_cell_params(fract_type *f)
{
  fprintf(stderr,"Input cell parameters a b c Alpha Beta Gamma \n");
  scanf("%lf %lf %lf %lf %lf %lf",&f->A,&f->B,&f->C,&f->Alpha,&f->Beta,&f->Gamma);
}

  

