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
  
  FILE : rdc3d.c
  AUTHOR(S) : Pat Walters
  DATE : 10-15-93
  PURPOSE :  Read a Chem3D file
  
  Modified 5-93 by PW  Allowed the use of single and multi-structure files
  Modified 12-93 by PW Allowed coordinates to be Cartesian or Fractional
  ******/

#include "bbltyp.h"

int read_mmads(FILE *file1, ums_type *mol)
{
  read_chem3d(file1,mol,"MMADS","MM2");
  print_ums(mol);
  return(TRUE);
}

int read_chem3d1(FILE *file1, ums_type *mol)
{
  read_chem3d(file1,mol,InputKeywords,"MM2");
  return(TRUE);
}

int read_chem3d2(FILE *file1, ums_type *mol)
{
  read_chem3d(file1,mol,InputKeywords,"C3D");
  return(TRUE);
}


int 
  read_chem3d(FILE *file1, ums_type *mol, char *keywords,char *type_key)
{
  char the_line[BUFF_SIZE];
  int i,k;
  char temp_type[5];
  char the_token[10];
  int tokens;
  int *label;
  char atomic_type[10];
  matrix_3x3 m;
  fract_type f;
  double exponent = 0.0;
  double divisor = 1.0;
  int has_fractional = FALSE, has_divisor = FALSE;
  int column;
  
  fgets(the_line,sizeof(the_line),file1);
  tokens = count_tokens(the_line," \t\n");

  if (EQ(keywords,"MMADS"))
  {
    for (i = 0; i < (int) strlen(the_line); i++)
      if ((!isdigit(the_line[i])) && (!isspace(the_line[i])))
      {
	the_line[i] = '\0';
	break;
      }
    Atoms = atoi(the_line);
  }
  
  switch(tokens)
  {
  case 7 :
    sscanf(the_line,"%d%lf%lf%lf%lf%lf%lf",
	   &Atoms,&f.Alpha,&f.Beta,&f.Gamma,&f.A,&f.B,&f.C);
    fill_orth_matrix(&f,&m);
    has_fractional = TRUE;
    break;
  case 8 :
    sscanf(the_line,"%d%lf%lf%lf%lf%lf%lf%lf",
	   &Atoms,&f.Alpha,&f.Beta,&f.Gamma,&f.A,&f.B,&f.C,&exponent);
    fill_orth_matrix(&f,&m);
    has_fractional = TRUE;
    has_divisor = TRUE;
    break;   
  default :
    sscanf(the_line,"%d",&Atoms);
    break;
  }

  divisor = pow(10.0,exponent);
  
  initialize_ums(&mol);

  ShowProgress(Atoms,"Reading Atoms");

  label = (int *)malloc((Atoms + 1) * sizeof(int));

  sscanf(the_line,"%d",&Atoms);

  column = locate_input_type(type_key);
  for (i = MIN_ATOM; i <= Atoms; i++)
  {

      UpdateProgress();
      fgets(the_line,sizeof(the_line),file1);
      sscanf(the_line,"%s%d%lf%lf%lf%s",
	     atomic_type,
	     &label[i],
	     &X(i),
	     &Y(i),
	     &Z(i),
	     temp_type);
      if (has_fractional) 
	fract_to_cart(&Point(i),&m);
      if (has_divisor)
      {
	X(i) = X(i)/divisor;
	Y(i) = Y(i)/divisor;
	Z(i) = Z(i)/divisor;
      }
      
      tokens = count_tokens(the_line,"\t\n ");
      for (k = 7; k <= tokens; k++)
      {
	strcpy(the_token,gettoken(the_line,"\t\n ",k));
	Connection(i,(k-7)) = atoi(the_token);
      }
      Atomic_number(i) = get_input_type(i,column,temp_type,Type(i),zero);
    }
  
  figure_valence(mol);
  xlate_c3d_labels(mol,label);
  build_connection_table(mol);
  assign_bond_order(mol);
  
  free(label);
  return(TRUE);
}


void xlate_c3d_labels(ums_type *mol, int *label)
{
  int i,j,k;
  
  for (i = 1; i <= Atoms; i++)
  {
    for (j = 0; j < Valence(i); j++)
      for (k = 1; k <= Atoms; k++)
      {
	if (label[k] == Connection(i,j))
	{
	  Connection(i,j) = k;
	  break;
	}
      }
  }
}
