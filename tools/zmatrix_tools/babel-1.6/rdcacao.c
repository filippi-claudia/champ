/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
---------------------------------------------------------------------------
FILE : read_cacao.c
AUTHOR(S) : Matt Stahl
DATE : 7-93
PURPOSE : routines to read the format used by the program Cacao
******/

#include "bbltyp.h"

#define DELIMS "\t\n ,"

int 
read_caccrt(FILE *file1, ums_type *mol)
{
  char cacao_line[BUFF_SIZE];
  int i;
  int result;
  fract_type f;
  matrix_3x3 m;
  
  fgets(cacao_line,sizeof(cacao_line), file1);
  fgets(cacao_line,sizeof(cacao_line), file1);

  sscanf(&cacao_line[0],"%d",&Atoms);
  
#ifdef MAC
ShowProgress(Atoms,"Reading Atoms");
#endif

  result = initialize_ums(&mol);
  
  while ( fgets(cacao_line,sizeof(cacao_line), file1) != NULL &&
         NOTEQn(cacao_line,"CELL",4) )
  {
  }
  
  sscanf(cacao_line,"%*s%*s%lf%lf%lf%lf%lf%lf",
	 &f.A,&f.B,&f.C,&f.Alpha,&f.Beta,&f.Gamma); 
  f.A = atof(gettoken(cacao_line,DELIMS,2));
  f.B = atof(gettoken(cacao_line,DELIMS,3));
  f.C = atof(gettoken(cacao_line,DELIMS,4));
  f.Alpha = atof(gettoken(cacao_line,DELIMS,5));
  f.Beta = atof(gettoken(cacao_line,DELIMS,6));
  f.Gamma = atof(gettoken(cacao_line,DELIMS,7));
  
  fill_orth_matrix(&f,&m);

  for (i = 1; i <= Atoms;i ++)
    {
#ifdef MAC
      UpdateProgress();
#endif
      fgets(cacao_line,sizeof(cacao_line), file1);
      sscanf(cacao_line,"%s",Type(i));
      X(i) = atof(gettoken(cacao_line,DELIMS,2));
      Y(i) = atof(gettoken(cacao_line,DELIMS,3));
      Z(i) = atof(gettoken(cacao_line,DELIMS,4));
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

    
	  
