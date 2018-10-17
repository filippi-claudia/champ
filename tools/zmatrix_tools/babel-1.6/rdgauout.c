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

FILE : rdgauout.c
AUTHOR(S) : Pat Walters
DATE : 1-94
PURPOSE : routines to read a Gaussian Output file
******/


#include "bbltyp.h"
#define DELIMS "\n\t ()"

int 
read_gau_out(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i;
  long pos = 0;
  int result;
  double the_energy = 0;
  int is_cartesian = FALSE;

  Atoms = 0;

  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (strstr(the_line,"GradGradGrad") && pos)
      break;

    if (strstr(the_line,"SCF DONE"))
    {
      the_energy = atof(gettoken(the_line,DELIMS,6));
      the_energy *= 627.5095;
    }
    
    if (strstr(the_line,"Z-MATRIX (ANGSTROMS AND DEGREES)") != NULL)
    {
      for (i = 0; i < 2; i++)
	fgets(the_line,sizeof(the_line), file1);
      
      pos = ftell(file1);
      fgets(the_line,sizeof(the_line), file1);
      while (!strstr(the_line,"-----------------------"))
      {
	Atoms++;
	fgets(the_line,sizeof(the_line), file1);
      }
    }
  }

  if (!file1)
    {
      Atoms = 0;
      return(FALSE);
    }
  
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  initialize_internal(&mol);
  Energy = the_energy;


  /* Need to check if the output coordiantes are Cartesian or
     Internal */

  fseek(file1,pos,0);
  fgets(the_line,sizeof(the_line), file1);
  if (count_tokens(the_line,DELIMS) == 7)
    is_cartesian = TRUE;
  fseek(file1,pos,0);
  
  /* Read the coordiantes */

  if (is_cartesian)
    read_gaussian_cartesian_atoms(file1,mol);
  else
    read_gaussian_internal_atoms(file1,mol);

  /* Make this into a UMS */
  
  if (Atoms > 0)
  {
    if (!is_cartesian)
      result = int_to_cart(mol);
    result = assign_radii(mol);
    result = assign_bonds(mol);
    result = assign_types(mol);
    result = build_connection_table(mol);
    assign_bond_order(mol);
  }
  
  pos = ftell(file1);
  while (fgets(the_line,sizeof(the_line), file1))
  {
    if (strstr(the_line,"Z-MATRIX (ANGSTROMS AND DEGREES)") != NULL)
    {
      while (fgets(the_line,sizeof(the_line), file1))
	{
	  if (strstr(the_line,"GradGradGrad") != NULL)
	    {
	      fseek(file1,pos,0);
	      return(TRUE);
	    }
	}
    }
    pos = ftell(file1);
  }
  
  return(TRUE);
}



void read_gaussian_cartesian_atoms(FILE *file1, ums_type *mol)
{
  int i;
  char the_line[BUFF_SIZE];
  
  for (i = 1; i <= Atoms; i++)
  {
    fgets(the_line,sizeof(the_line), file1);
    sscanf(the_line,"%*s%*s%s%*s%lf%lf%lf",Type(i),&X(i),&Y(i),&Z(i));
  }
}

void read_gaussian_internal_atoms(FILE *file1, ums_type *mol)
{
  int i;
  char the_line[BUFF_SIZE];
  
  i = 0;
  while (i < Atoms)
  {
    UpdateProgress();
    fgets(the_line,sizeof(the_line), file1);
    i++;
/*
    printf("%d -> %s\n",i,the_line); 
*/

    switch(i)
    {
    case 1 :
      strcpy(Type(1),gettoken(&the_line[9],DELIMS,1));
      break;
      
    case 2 :
      strcpy(Type(2),gettoken(&the_line[9],DELIMS,1));
      NA(2) = atoi(gettoken(&the_line[13],DELIMS,1));
      R(2) = atof(gettoken(&the_line[13],DELIMS,2));
      break;
      
    case 3 :
      strcpy(Type(3),gettoken(&the_line[9],DELIMS,1));
      NA(3) = atoi(gettoken(&the_line[13],DELIMS,1));
      R(3) = atof(gettoken(&the_line[13],DELIMS,2));
      NB(3) = atoi(gettoken(&the_line[13],DELIMS,4));
      W(3) = atof(gettoken(&the_line[13],DELIMS,5));
      break;
      
    default :
      strcpy(Type(i),gettoken(&the_line[9],DELIMS,1));
      NA(i) = atoi(gettoken(&the_line[13],DELIMS,1));
      R(i) = atof(gettoken(&the_line[13],DELIMS,2));
      NB(i) = atoi(gettoken(&the_line[13],DELIMS,4));
      W(i) = atof(gettoken(&the_line[13],DELIMS,5));
      NC(i) = atoi(gettoken(&the_line[13],DELIMS,7));
      T(i) = atof(gettoken(&the_line[13],DELIMS,8));
      break;
    }
    clean_atom_type(Type(i));
  }
}




   
    
    
    
	  


