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
FILE : rdpsvbin.c
AUTHOR(S) : Pat Walters
DATE : 12-93
PURPOSE : routines to read a Gaussian Z-matrix
******/

#include "bbltyp.h"
#define PSGVB_DELIMS "\t\n, =@#*"

#define Sym(x) the_syms[x].sym
#define Val(x) the_syms[x].val

int
read_psgvb_input(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i = 0, j;
  int sym_count = 0;
  int result;
  zsymbol *the_syms;
  char lab1[10],lab2[10],lab3[10];
  char sym1[30],sym2[30],sym3[30];
  int blank_lines = 0;
  int stop = FALSE;
  long int atom_start = -1;
  long int variable_start = -1;
  char first_token[BUFF_SIZE];
  int tokens;
  double factor;
  int itok;
  int zmat_style = 0;
#ifdef DEBUG
  int iline = 0;
#endif

  Atoms = 0;
  /* Next section */
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
  {
#ifdef DEBUG
    printf("%d -> %s\n",++iline,the_line);
#endif
    if (count_tokens(the_line,PSGVB_DELIMS) == 0)
      continue;
    strcpy(first_token,gettoken(the_line,PSGVB_DELIMS,1));
    lowercase(first_token);
    if (EQ(first_token,"$zmat") || EQ(first_token,"&zmat"))
      {
        atom_start = ftell(file1);
#ifdef DEBUG
        printf("(atom_start)\n");
#endif
        i = 0;
        while ( fgets(the_line,sizeof(the_line), file1) != NULL)
          {
#ifdef DEBUG
            printf("%d (Atom %d) -> %s\n",++iline,i+1,the_line);
#endif
            if (the_line[0] == '$' || the_line[0] == '&')
              break;

            i++;
            Atoms++;
            if (i == 1 && count_tokens(the_line,PSGVB_DELIMS) == 1)
              zmat_style=1;
          }
      }
    if (EQ(first_token,"$zvar") || EQ(first_token,"&zvar"))
      {
        variable_start = ftell(file1);
#ifdef DEBUG
        printf("(variable_start)\n");
#endif
      }
  }

  initialize_ums(&mol);
  initialize_internal(&mol);

  /* Read zvars */
  the_syms = (zsymbol *)malloc(3 * Atoms * sizeof(zsymbol));
  sym_count = 0;
  fseek(file1,variable_start,0);
  while (fgets(the_line,sizeof(the_line), file1) != NULL)
    {
#ifdef DEBUG
      printf("%d -> %s",sym_count,the_line);
#endif
      if (the_line[0] == '$' || the_line[0] == '&')
        break;

      tokens = count_tokens(the_line,PSGVB_DELIMS);
      itok = 0;
      while ( tokens >= 2 )
        {
          strcpy(Sym(sym_count),gettoken(the_line,PSGVB_DELIMS,++itok));
          Val(sym_count) = (double)atof(gettoken(the_line,PSGVB_DELIMS,++itok));
#ifdef DEBUG
          printf("  %d   %s = %f\n",sym_count,Sym(sym_count),Val(sym_count));
#endif
          tokens -= 2;
          sym_count++;
        }
    }

  /* Read zmat */
  ShowProgress(Atoms,"Reading Atoms");
  fseek(file1,atom_start,0);
  for (j = 1; j <= Atoms; j++)
    {
      UpdateProgress();
      fgets(the_line,sizeof(the_line), file1);
      tokens = count_tokens(the_line,PSGVB_DELIMS);
      if (zmat_style == 1)
        {
          switch(j)
            {
            case 1 :
              strcpy(Type(j),gettoken(the_line,PSGVB_DELIMS,1));
              break;
            case 2 :
              strcpy(Type(j),gettoken(the_line,PSGVB_DELIMS,1));
              strcpy(lab1,gettoken(the_line,PSGVB_DELIMS,2));
              strcpy(sym1,gettoken(the_line,PSGVB_DELIMS,3));

              NA(j) = xlate_label(mol,lab1,j);
              R(j) = xlate_ps_symbol(the_syms,sym1,sym_count);
              break;
            case 3 :
              sscanf(the_line,"%s,%s,%s,%s,%s",Type(j),lab1,sym1,lab2,sym2);
              strcpy(Type(j),gettoken(the_line,PSGVB_DELIMS,1));
              strcpy(lab1,gettoken(the_line,PSGVB_DELIMS,2));
              strcpy(sym1,gettoken(the_line,PSGVB_DELIMS,3));
              strcpy(lab2,gettoken(the_line,PSGVB_DELIMS,4));
              strcpy(sym2,gettoken(the_line,PSGVB_DELIMS,5));

              NA(j) = xlate_label(mol,lab1,j);
              NB(j) = xlate_label(mol,lab2,j);
              R(j) = xlate_ps_symbol(the_syms,sym1,sym_count);
              W(j) = xlate_ps_symbol(the_syms,sym2,sym_count);
              break;
            default :
              strcpy(Type(j),gettoken(the_line,PSGVB_DELIMS,1));
              strcpy(lab1,gettoken(the_line,PSGVB_DELIMS,2));
              strcpy(sym1,gettoken(the_line,PSGVB_DELIMS,3));
              strcpy(lab2,gettoken(the_line,PSGVB_DELIMS,4));
              strcpy(sym2,gettoken(the_line,PSGVB_DELIMS,5));
              strcpy(lab3,gettoken(the_line,PSGVB_DELIMS,6));
              strcpy(sym3,gettoken(the_line,PSGVB_DELIMS,7));

              NA(j) = xlate_label(mol,lab1,j);
              NB(j) = xlate_label(mol,lab2,j);
              NC(j) = xlate_label(mol,lab3,j);
              R(j) = xlate_ps_symbol(the_syms,sym1,sym_count);
              W(j) = xlate_ps_symbol(the_syms,sym2,sym_count);
              T(j) = xlate_ps_symbol(the_syms,sym3,sym_count);
              break;
            }
        }
      else
        {
          strcpy(Type(j),gettoken(the_line,PSGVB_DELIMS,1));
          strcpy(sym1,gettoken(the_line,PSGVB_DELIMS,2));
          strcpy(sym2,gettoken(the_line,PSGVB_DELIMS,3));
          strcpy(sym3,gettoken(the_line,PSGVB_DELIMS,4));

          X(j) = xlate_ps_symbol(the_syms,sym1,sym_count);
          Y(j) = xlate_ps_symbol(the_syms,sym2,sym_count);
          Z(j) = xlate_ps_symbol(the_syms,sym3,sym_count);
        }
    }
  for (i = 1; i <= Atoms; i++)
    clean_atom_type(Type(i));

  if (Atoms > 0)
    {
      if (zmat_style == 1)
        {
          result = int_to_cart(mol);
        }
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


double xlate_ps_symbol(zsymbol *the_syms, char *sym, int max)
{
  int i;
  int is_neg = FALSE;
  char the_str[100];
  char numeric[] = "0123456789.+";

  if (strchr(numeric,sym[0]) != NULL ||
      (sym[0] == '-' && strchr(numeric,sym[1]) != NULL))
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

