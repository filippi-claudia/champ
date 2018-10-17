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

FILE : rdshelx.c
AUTHOR(S) : Pat Walters
DATE : 12 -93
PURPOSE : routines to read a Shelx file
******/

#include "bbltyp.h"

int 
read_shelx(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int found = FALSE;
  fract_type f;
  int i = 0;
  int result;
  matrix_3x3 m;

  Atoms = 0;
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (count_tokens(the_line,"\n\t ") > 0)
      if EQ(gettoken(the_line,"\n\t ",1),"CELL")
      {
	found = TRUE;
	sscanf(the_line,"%*s%*s%lf%lf%lf%lf%lf%lf",
	       &f.A,&f.B,&f.C,&f.Alpha,&f.Beta,&f.Gamma);
	fill_orth_matrix(&f,&m);
      }
      else
	if (found)
	{
	  if (is_good_shelx_line(the_line))
	    Atoms++;
	}      
  }
  ShowProgress(Atoms,"Reading Atoms");
  result = initialize_ums(&mol);
  rewind(file1);
  i = 0;
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
    if (is_good_shelx_line(the_line))
    {
      UpdateProgress();
      i++;
      sscanf(the_line,"%s %*s %lf %lf %lf",
	     Type(i),
	     &X(i),
	     &Y(i),
	     &Z(i));
      check_shelx_coords(&Point(i)); 
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

int count_shelx_atoms(char *the_line)
{
  int atom_count = 0;
  int i;
  int tokens;
  char the_token[20];
  
  tokens = count_tokens(the_line,"\t\n ");
  for (i = 1; i <= tokens; i++)
  {
    strcpy(the_token,gettoken(the_line,"\t\n ",i));
    if (isdigit(the_token[0]))
      atom_count += atoi(the_token);
  }
  return(atom_count);
}

	   
void check_shelx_coords(coord_type *p)
{
  if (p->x > 10.0)
    p->x -= 10.0;
  if (p->y > 10.0)
    p->y -= 10.0;
  if (p->z > 10.0)
    p->z -= 10.0;
}


int is_good_shelx_line(char *the_line)
{
  char first[BUFF_SIZE];
  int possible = FALSE;
  int has_digit = FALSE;
  int i;

  if ((strchr(the_line,'(')) || (strchr(the_line,')')))
    return(FALSE);

  if ((count_tokens(the_line,"\n\t ") >= 4) && (isalpha(the_line[0])))
  {
    strcpy(first,gettoken(the_line," \n\t",1));
    if (isdigit(first[0]))
      return(FALSE);
    for (i = 0; i < (int) strlen(first); i++)
    {
      if (isdigit(first[i]))
      {
	has_digit = TRUE;
	break;
      }
    }
    if ((!has_digit) && (strlen(first) > 2))
      return(FALSE);
    clean_atom_type(first);
    if (is_element(first))
    {
      return(TRUE);
    }
  }
  return(FALSE);
}

