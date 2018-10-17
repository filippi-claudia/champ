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

FILE : typbybo.c
AUTHOR(S) : Pat Walters
DATE : 2-96
PURPOSE : assign types according to bond order

******/

#include "bbltyp.h"

#undef DEBUG

void assign_type_by_bo(ums_type *mol)
{
  int i,j, conn;
  int max_order;
  int free_ox;
  int is_isocyanate;

  for (i = 1;i <= Atoms;i++)
    {
      Atomic_number(i) = get_atomic_number(Type(i));

      max_order = 0;
      for (j = 0;j < Valence(i);j++)
	if (BO(i,j) > max_order)
	  max_order = BO(i,j);
      
      switch (Atomic_number(i))
	{
	case 6:   /* C */
	case 7:   /* N */
	case 8:   /* O */
	case 15:  /* P */
	case 16:  /* S */
	
	  switch (max_order)
	    {
	    case 1:
	      strcat(Type(i),"3");
	      break;
	    case 5:
	    case 2:
	      strcat(Type(i),"2");
	      break;
	    case 3:
	      strcat(Type(i),"1");
	      break;	
	    }
	}

      if ((Atomic_number(i) == 7) && ((Valence(i) == 4) || (Valence(i) == 0)))
      {
	if (count_free_ox(mol,i) >= 1)
	  strcpy(Type(i),"Nox");
	else
	  strcpy(Type(i),"N3+");
      }

      if (EQ(Type(i),"O1"))
	strcpy(Type(i),"O2");

      if (EQ(Type(i),"N2"))
	strcpy(Type(i),"Npl");
    }

  type_hydrogens(mol);

  /* assign Cac Pac Sac Nox and O- to attached oxygens*/

  for (i = 1;i <= Atoms;i++)
    {
      free_ox = count_free_ox(mol,i);
#ifdef DEBUG
      printf("i = %d type = %s valence = %d free_ox = %d atomic_number = %d heavy = %d \n",
	     i,Type(i),Valence(i),free_ox,Atomic_number(i),count_heavy_atoms(mol,i));
#endif
      switch (Atomic_number(i))
	{
	case 6:   /* C */
	  if (EQ(Type(i),"C2") && (free_ox >= 2))
	    {
	      strcpy(Type(i),"Cac");
	      type_attached_oxygens(mol,i);
	    }
	  
	  /* Type isocyanates as C1 */

	  if (EQ(Type(i),"C2") && (Valence(i) == 2))
	  {
	    is_isocyanate = TRUE;
	    for (j = 0; j < Valence(i); j++)
	    {
	      if (BO(i,j) != 2)
		is_isocyanate = FALSE;
	    }
	    if (is_isocyanate)
	      strcpy(Type(i),"C1");	    
	  }
	  break;
	case 7:  /* N */
	  if (EQ(Type(i),"N3"))
	  {
	    /* Nitrogens next to a unsaturated atom or a planar nitrogen (hydrazines) 
	       are assigned as sp2 */

	    for (j = 0; j < Valence(i); j++)
	    {
	      conn = Connection(i,j);
	      if (IsUnsatType(Type(conn)) || EQ(Type(conn),"Npl") || EQ(Type(conn),"N2"))
	      {
		strcpy(Type(i),"Npl");
		break;
	      }
	    }
	  }

	  if (EQ(Type(i),"Npl") && (free_ox >= 2))
	    {
	      strcpy(Type(i),"Ntr");
	      type_attached_oxygens(mol,i);
	    }

      	  if (EQ(Type(i),"Npl"))
	  {
	    for (j = 0; j < Valence(i); j++)
	      if (BO(i,j) == 2)
		strcpy(Type(i),"N2");
	  }	    

	  if (EQ(Type(i),"Nox") && (free_ox >= 2))
	    {
	      type_attached_oxygens(mol,i);
	    }
	  break;
	case 15:  /* P */
	  if (Valence(i) == 4)
	  {
	    if (free_ox >= 2)
	    {
	      strcpy(Type(i),"Pac");
	      type_attached_oxygens(mol,i);
	    }
	    else if (free_ox == 1)
	    {
	      strcpy(Type(i),"Pox");
	      type_attached_oxygens(mol,i);
	    }
	    else
	      strcpy(Type(i),"P3+");
	  }
	  else
	    strcpy(Type(i),"P3");
	  break;
	case 16:  /* S */
	  if ((free_ox >= 3) || ((Valence(i) == 3) && (free_ox == 2)))
	    {
	      strcpy(Type(i),"Sac");
	      type_attached_oxygens(mol,i);
	    }
	  else
	    if (free_ox == 2)
	    {
	      strcpy(Type(i),"So2");
	      type_attached_oxygens(mol,i);
	    }
	    else
	      if (free_ox >= 1)
	      {
		strcpy(Type(i),"Sox");
		type_attached_oxygens(mol,i);
	      }
	  break;
	}
    }
  phase6(mol);
  check_for_amides(mol);
}


void type_attached_oxygens(ums_type *mol, int atm)
{
  int i, conn;
  int is_O_minus = FALSE;
  

  if (EQ(Type(atm),"Cac") || EQ(Type(atm),"Sac") || EQ(Type(atm),"Pac"))
  {
    is_O_minus = TRUE;
  }

  for (i = 0; i < Valence(atm); i++)
    {
      conn = Connection(atm,i);
#ifdef DEBUG
      printf("atm = %d i = %d conn = %d Atomic number = %d Valence = %d \n",
	     atm,i,conn,Atomic_number(conn),Valence(conn));
#endif
      if ((Atomic_number(conn) == 8) && (count_heavy_atoms(mol,conn) == 1))
      {
	if (is_O_minus)
	  strcpy(Type(conn),"O-");
	else
	  strcpy(Type(conn),"O2");
      }
    }
}


