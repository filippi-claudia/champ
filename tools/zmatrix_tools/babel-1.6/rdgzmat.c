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
FILE : rdgzmat.c
AUTHOR(S) : Pat Walters
DATE : 12-93
PURPOSE : routines to read a Gaussian Z-matrix
******/

#include "bbltyp.h"
#define GAU_DELIMS "\t\n, "

int 
read_gau_zmatrix(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i = 0, j;
  int sym_count = 0;
  int result;
  zsymbol *the_syms;
  char lab1[10],lab2[10],lab3[10];
  char sym1[10],sym2[10],sym3[10];
  int blank_lines = 0;
  int stop = FALSE;
  long int atom_start,variable_start;
  char first_token[BUFF_SIZE];
  int tokens;
  double factor;

  Atoms = 0;
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (the_line[0] == '#')
      break;
  }
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (is_blank_line(the_line))
      blank_lines++;
    if (blank_lines == 2)
      break;
  }
  toss(file1,1);
  atom_start = ftell(file1);
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
    {
      if (is_blank_line(the_line))
	stop = TRUE;
      else
      {
	strcpy(first_token,gettoken(the_line,", \t\n",1));
	uppercase(first_token);
	if (EQn(first_token,"VARIABLE",8))
	  stop = TRUE;
      }
      if (stop == TRUE)
	break;
      else
	Atoms++;
    }
  initialize_ums(&mol);
  initialize_internal(&mol);
  variable_start = ftell(file1);
  i = 0;
  while ((fgets(the_line,sizeof(the_line), file1) != NULL))
  {
    i++;
  }
  the_syms = (zsymbol *)malloc(3 * Atoms * sizeof(zsymbol));
  fseek(file1,variable_start,0);
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
    if (count_tokens(the_line," =\t\n") == 2)
    {
      strcpy(the_syms[sym_count].sym,gettoken(the_line," =\t\n",1));
      the_syms[sym_count].val = (double)atof(gettoken(the_line," =\t\n",2)); 
      sym_count++;
    }
  }
  ShowProgress(Atoms,"Reading Atoms");  
  fseek(file1,atom_start,0);
  for (j = 1; j <= Atoms; j++)
  {
    factor = 0.0;
    UpdateProgress();
    fgets(the_line,sizeof(the_line), file1);
    tokens = count_tokens(the_line,GAU_DELIMS);
    switch(j)
      {
      case 1 :
	strcpy(Type(j),gettoken(the_line,GAU_DELIMS,1));
	break;
      case 2 :
	strcpy(Type(j),gettoken(the_line,GAU_DELIMS,1));
	strcpy(lab1,gettoken(the_line,GAU_DELIMS,2));
	strcpy(sym1,gettoken(the_line,GAU_DELIMS,3));

	mol->internal[j].na = xlate_label(mol,lab1,j);
	mol->internal[j].r = xlate_symbol(the_syms,sym1,sym_count);
	break;
      case 3 :
	sscanf(the_line,"%s,%s,%s,%s,%s",Type(j),lab1,sym1,lab2,sym2);
	strcpy(Type(j),gettoken(the_line,GAU_DELIMS,1));
	strcpy(lab1,gettoken(the_line,GAU_DELIMS,2));
	strcpy(sym1,gettoken(the_line,GAU_DELIMS,3));
	strcpy(lab2,gettoken(the_line,GAU_DELIMS,4));
	strcpy(sym2,gettoken(the_line,GAU_DELIMS,5));

	mol->internal[j].na = xlate_label(mol,lab1,j);
	mol->internal[j].nb = xlate_label(mol,lab2,j);
	mol->internal[j].r = xlate_symbol(the_syms,sym1,sym_count);
	mol->internal[j].w = xlate_symbol(the_syms,sym2,sym_count);
	break;
      default :
	strcpy(Type(j),gettoken(the_line,GAU_DELIMS,1));
	strcpy(lab1,gettoken(the_line,GAU_DELIMS,2));
	strcpy(sym1,gettoken(the_line,GAU_DELIMS,3));
	strcpy(lab2,gettoken(the_line,GAU_DELIMS,4));
	strcpy(sym2,gettoken(the_line,GAU_DELIMS,5));
	strcpy(lab3,gettoken(the_line,GAU_DELIMS,6));
	strcpy(sym3,gettoken(the_line,GAU_DELIMS,7));
	if (tokens == 8)
	  factor = atof(gettoken(the_line,GAU_DELIMS,8));

	mol->internal[j].na = xlate_label(mol,lab1,j);
	mol->internal[j].nb = xlate_label(mol,lab2,j);
	mol->internal[j].nc = xlate_label(mol,lab3,j);
	mol->internal[j].r = xlate_symbol(the_syms,sym1,sym_count);
	mol->internal[j].w = xlate_symbol(the_syms,sym2,sym_count);
	mol->internal[j].t = xlate_symbol(the_syms,sym3,sym_count);
	if (factor != 0.0)
	  mol->internal[j].t *= factor;
	break;
      }
  }
  for (i = 1; i <= Atoms; i++)
    clean_atom_type(Type(i));

  if (Atoms > 0)
  {
    result = int_to_cart(mol);
    result = assign_radii(mol);
    result = assign_bonds(mol);
    result = assign_types(mol);
    result = build_connection_table(mol);
    assign_bond_order(mol);
  }
  free(the_syms);
  read_to_eof(file1);
  return(TRUE);
}



int xlate_label(ums_type *mol, char *sym, int max)
{
  int i;
  char the_str[100];
  
  if (isdigit(sym[0]))
    return(atoi(sym));
  for (i = 1; i < max; i++)
    if (EQ(Type(i),sym))
      return(i);
  sprintf(the_str,"cant translate label %s\n",sym);
  show_warning(the_str);
  return(0);
}
  
double xlate_symbol(zsymbol *the_syms, char *sym, int max)
{
  int i;
  int is_neg = FALSE;
  char the_str[100];

  if ((sym[0] == '-') ||  (isdigit(sym[0])))
    return(atof(sym));
  
  if (sym[0] == '-')
  {
    is_neg = TRUE;
    sym++;
  }

  for (i = 0; i < max; i++)
    if (EQ(the_syms[i].sym,sym))
    {
      if (is_neg == TRUE)
	return(the_syms[i].val * (-1));
      else
	return(the_syms[i].val);
    }
  sprintf(the_str,"cant translate symbol %s\n",sym);
  show_warning(the_str);
  return(0.0);
}





   
    
    
    
	  


