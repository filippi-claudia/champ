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
FILE : addh2.c
AUTHOR(S) : Pat Walters
DATE : 2-6-96

New hydrogen addition code.  Uses bond orders to figure out how many
hydrogens to add.
******/

#undef DEBUG

#include "bbltyp.h"

int count_missing_bo_hydrogens(ums_type *mol)
{
  int missing = 0;
  int i;
  char temp_type[5];
  int to_add, type_valence, attached;
  int result;
  int done;

  for (i = 1; i <= Atoms; i++)
    {
      type_valence = 0;

      switch(Atomic_number(i))
	{
	case 6 :
	  type_valence = 4;
	  break;
	case 7 :
	  if ((EQ(Type(i),"N2")) && (Valence(i) == 1))
	    type_valence = 2;
	  else
	    if (EQ(Type(i),"N3+"))
	      type_valence = 4;
	    else
	      type_valence = 3;
	  break;
	case 8 :
	  if ((EQ(Type(i),"O-")) || (EQ(Type(i),"O2")))
	    type_valence = 1;
	  else
	    type_valence = 2;
	  break;
	}
  
      attached = count_attached_bonds(mol,i);
      
      Redo(i) = 0;
      to_add = 0;
      if ((Valence(i) < type_valence) && (Valence(i) > 0))
	{
	  to_add = type_valence - attached;
	  if (to_add > 0)
	    {
	      Redo(i) = to_add;
	      missing += to_add;
	    }
	}
  
#ifdef DEBUG
      printf("num = %d type = %s atomic_number = %d attached = %d type_valence = %d to_add = %d\n",
	     i,Type(i),Atomic_number(i),attached,type_valence,to_add);
#endif
    }
  return(missing);
}


void place_hydrogens2(ums_type *mol, int old_count, int num_H_to_add)
{
  char type_name[5]; /* string to hold hybridization */
  int hyb;           /* hybridization 0,1,2,3 */
  int code;          /* atom type code 10 * atomic number + hybridization */
  int to_add;        /* number of hydrogens to add to each atom */
  int h_count;       /* number of hydrogens added to the molecule */
  char temp_title[BUFF_SIZE];
  int i;

  Atoms += num_H_to_add;
  h_count = old_count + 1;
  strcpy(temp_title,Title);
  reinitialize_ums(&mol);
  zero_out_ums(mol,old_count + 1);
  strcpy(Title,temp_title);

  for (i = 1; i <= old_count; i++)
    {
      get_output_type(i,"HYB",Type(i),type_name,zero);
      hyb = atoi(type_name);
      code = Atomic_number(i) * 10 + hyb;
      to_add = Redo(i);

#ifdef DEBUG
      printf("i = %d code = %d to_add = %d h_count = %d\n",
	     i,code,to_add,h_count);
#endif

      switch(code)
	{
	case 63 :   /* sp3 carbon */
	  switch(to_add)
	    {
	    case 1 :
	      add_tertiary_hydrogen(mol,i,h_count,SP3_C_H_DIST);
	      h_count += 1;
	      break;
	    case 2 :
	      add_methylene_hydrogens(mol,i,h_count,SP3_C_H_DIST);
	      h_count += 2;
	      break;
	    case 3 :
	      add_methyl_hydrogen(mol,i,h_count,SP3_C_H_DIST);
	      h_count += 1;
	      add_methylene_hydrogens(mol,i,h_count,SP3_C_H_DIST);
	      h_count += 2;
	      break;
	    }
	  break;
	case 73 : /* sp3 nitrogen */
	  switch(to_add)
	    {
	    case 1 :
	      if (EQ(Type(i),"N3+"))
		add_tertiary_hydrogen(mol,i,h_count,SP3_N_H_DIST);
	      else
		add_sp3_N_hydrogen(mol,i,h_count,SP3_N_H_DIST);
	      h_count += 1;
	      break;
	    case 2 :
	      add_methylene_hydrogens(mol,i,h_count,SP3_N_H_DIST);
	      h_count += 2;
	      break;
	    case 3 :
	      add_methyl_hydrogen(mol,i,h_count,SP3_N_H_DIST);
	      h_count += 1;
	      add_methylene_hydrogens(mol,i,h_count,SP3_N_H_DIST);
	      h_count += 2;
	      break;
	    }
	  break;
	case 62 :  /* sp2 carbon */
	  switch(to_add)
	    {
	    case 1 :
	      add_sp2_hydrogen(mol,i,h_count,SP2_C_H_DIST);
	      h_count += 1;
	      break;
	    }
	  break;
	case 72 : /* sp2 nitrogen */
	  switch(to_add)
	    {
	    case 1 :
	      add_sp2_hydrogen(mol,i,h_count,SP2_N_H_DIST);
	      h_count += 1;
	      break;
	    }
	  break;
	case 61 : /* sp carbon */
	  switch(to_add)
	    {
	    case 1 :
	      add_sp_hydrogen(mol,i,h_count,SP_C_H_DIST);
	      h_count++;
	      break;
	    }
	  break;
	case 83 : /* sp3 oxygen */
	  switch(to_add)
	    {
	    case 1:
	      add_methyl_hydrogen(mol,i,h_count,SP3_O_H_DIST);
	      h_count++;
	      break;
	    }
	  break;
	}
    }

  /* save vinyl and amide protons for last, this way we know where to put them */

  for (i = 1; i <= old_count; i++)
    {
      get_output_type(i,"HYB",Type(i),type_name,zero);
      hyb = atoi(type_name);
      code = Atomic_number(i) * 10 + hyb;
      to_add = Redo(i);
      
      switch (code)
	{
	case 62 : /* sp2 carbon */
	  switch(to_add)
	    {
	    case 2 :
	      add_vinyl_hydrogens(mol,i,h_count,SP2_C_H_DIST);
	      h_count += 2;
	      break;
	    }
	  break; 
	case 72 : /* sp2 nitrogens */
	  switch(to_add)
	    {
	    case 2 :
	      add_vinyl_hydrogens(mol,i,h_count,SP2_N_H_DIST);
	      h_count += 2;
	      break;
	    }
	}
    }
  
  Atoms = h_count - 1;
  build_connection_table(mol);
}



  
    

int count_attached_bonds(ums_type *mol, int atm)
{
  int i;
  float bonds = 0.0;
  
  for (i = 0; i < Valence(atm); i++)
  {
    if (BO(atm,i) == 5)
      bonds += 1.5;
    else
      bonds += (float)BO(atm,i);
  }
  return((int)bonds);
}





