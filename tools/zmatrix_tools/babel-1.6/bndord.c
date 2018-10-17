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

FILE : bndord.c
AUTHOR(S) : Pat Walters
DATE : 2-10-93
PURPOSE : Assign bond orders based on atom types, clean up conjugated pi systems
******/

#include "bbltyp.h"

static warning wstr;

void
assign_bond_order(ums_type *mol)
{
  if (Atoms > 200)
    assign_bond_order1(mol);
  else
    assign_bond_order2(mol);
}

void
assign_bond_order1(ums_type *mol)
{
  int i, result;
  int sum_code;
  int column;
  char temp_type[10];
  int hyb, hyb_val[4] = {0,3,2,1};

  for (i = 1; i <= Atoms; i++)
  {
    get_output_type(i,"HYB",Type(i),temp_type,dummy);
    hyb = atoi(temp_type);
    Redo(i) = hyb_val[hyb];
  }

  for (i = 0; i < Bonds; i++)
  {
    sum_code = Redo(Start(i)) + Redo(End(i));
    switch(sum_code)
    {
    case 6 :
      Bond_order(i) = 3;
      break;
    case 4 :
      Bond_order(i) = 2;
      break;
    default :
      Bond_order(i) = 1;
    }
    if (is_carboxyl(mol,i))
      Bond_order(i) = 2;
    if (Bond_order(i) < 1 || Bond_order(i) > 3)
    {
      sprintf(wstr,"Bond %d - atoms %d - %d is wierd - Bond order is %d\n",
	     i,Start(i),End(i),Bond_order(i));
    show_warning(wstr);
    }
/*    printf("%s %s i = %d sum_code = %d bo = %d \n",Type(Start(i)),Type(End(i)),i,sum_code,Bond_order(i));*/
  }
/*  result = check_for_overflow(mol);*/
  result = check_for_conjugation(mol);
  dissect_connection_table(mol);
}

int is_carboxyl(ums_type *mol, int the_bond)
{
  int c_end = 0;
  int o_end = 0;
  int check = FALSE;
  
  if (EQ(Type(Start(the_bond)),"Cac") && (Type(End(the_bond))[0] == 'O'))
  {
    c_end = Start(the_bond);
    o_end = End(the_bond);
    check = TRUE;
  }

  if (EQ(Type(End(the_bond)),"Cac") && (Type(Start(the_bond))[0] == 'O'))
  {
    check = TRUE;
    c_end = End(the_bond);
    o_end = Start(the_bond);
  }

  if ((check) && (Valence(o_end) == 1))
    return(TRUE);
  else
    return(FALSE);
}

int check_for_overflow(ums_type *mol)
{
  int i;
  
  for (i = 0; i < Bonds; i++)
  {
    if ((Valence(Start(i)) == Max_bonds(Start(i)) ||
	 (Valence(End(i)) == Max_bonds(End(i)))))
      Bond_order(i) = 1;
  }
  return(TRUE);
}

int 
check_for_conjugation(ums_type *mol)
{
  int i,j;
  
  for (i = 0; i < Bonds; i++)
  {
    for (j = 0; j < Bonds; j++)
      if ((i != j) && (Bond_order(i) > 1))
	if (atom_in_common(mol,i,j) && (Bond_order(j) > 1))
	{
	  if ((Valence(Start(j)) > 1 && (Valence(End(j)) > 1)))
	  Bond_order(j) = 1;
	}
  }
  return(TRUE);
}
/******************************
atom_in_common
returns true if bond1 and bond2 have
an atom in common
*******************************/

int 
atom_in_common(ums_type *mol, int bond1, int bond2)
{
  if (Start(bond1) == Start(bond2)) return TRUE;
  if (End(bond1) == End(bond2)) return TRUE;
  if (Start(bond1) == End(bond2)) return TRUE;
  if (End(bond1) == Start(bond2)) return TRUE;
  return FALSE;
}


int assign_bond_code(char *input_type)
{
  static struct std_types
  {
    char *internal;
    int  output;
  }  type_tab[] =
  {
    {"C3",1},
    {"C2",2},
    {"C1",3},
    {"Cac",2},
    {"N3+",1},
    {"N3",1},
    {"Nam",2},
    {"Npl",2},
    {"N1",3},
    {"Nox",1},
    {"Ntr",2},
    {"Ng+",1},
    {"NC",2},
    {"O3",1},
    {"O2",2},
    {"O-",2},
    {"S3+",1},
    {"S3",1},
    {"S2",2},
    {"Sac",2},
    {"Sox",2},
    {"S",1},
    {"Bac",1},
    {"Box",1},
    {"B",1},
    {"Pac",2},
    {"Pox",2},
    {"P3+",1},
    {"P",1},
    {"HC",1},
    {"H",1},
    {"DC",1},
    {"D",1},
    {"F",1},
    {"Cl",1},
    {"Br",1},
    {"I",1},
    {"Ge",1},
    {"Sn",1},
    {"Pb",1},
    {"Se",1},
    {"Te",1},
    {"O",1},
    {"DONE",0},
    };
  int count = 0;

  while (strcmp(type_tab[count].internal,"DONE") != 0)
  {
    if (EQ(input_type,type_tab[count].internal))
      return(type_tab[count].output);
    count ++;
  }
  return(0);
}






