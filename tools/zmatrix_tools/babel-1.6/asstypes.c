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
FILE : asstypes.c 
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : Assign specific atom types (i.e. hybridization) to the atoms in the UMS

The code here is mine, but the ideas are
not. The majority of this program is a 
translation of a fortran program by 
Elaine Meng - UC San Francisco.
For information on the algorithms used
here see
E. Meng and R. Lewis, J. Comp. Chem., 12,
pp 891-898 (1991)
******/

#include "bbltyp.h"

#undef DEBUG



void type_attached_oxygens(ums_type *mol, int atm);

void check_atomic_numbers(ums_type *mol)
{
  int i, column;
  char temp_type[5];
  
  if (Atomic_number(1) == 0)
  {
    column = locate_input_type("INT");  
    for (i = 1; i <= Atoms; i++)
	Atomic_number(i) = get_input_type(i,column,Type(i),temp_type,dummy);
  }
}


int assign_types(ums_type *mol)
{
  int i, result;

  for (i = 1; i <= Atoms; i++)
    Atomic_number(i) = get_atomic_number(Type(i));

  tag_organics(mol);
  result = phase1(mol);
  result = phase2(mol);
  result = phase3(mol); 
  result = phase4(mol);
  result = phase5(mol); 
  result = phase6(mol);
  
  check_for_amides(mol);
  /*
    find_aromatic_atoms(mol);
    */  
  return(1);
}


void tag_organics(ums_type *mol)
{
  int i;
  
  for (i = 1; i <= Atoms; i++)
  {
    if (EQ(Type(i),"C") ||
	EQ(Type(i),"H") ||
	EQ(Type(i),"O") ||
	EQ(Type(i),"N") ||
	EQ(Type(i),"S") ||
	EQ(Type(i), "P"))
      Organic(i) = TRUE;
    else
      Organic(i) = FALSE;
  }
}

/********************************************************************
phase1 - type all hydrogens and deuterium by whether they are
attached to carbon or not. Calculate the number of heavy atoms
bonded to each atom by subtracting the number of hydrogens 
attached from the valenece
*********************************************************************/

int 
phase1(ums_type *mol)
{
  int result;
  result = type_hydrogens(mol);
  return(1);
}


int 
phase2(ums_type *mol)
{
  int result;

  result = valence_four(mol);
  result = valence_three(mol);
  result = valence_two(mol);
  return(1);
}

int 
phase3(ums_type *mol)
{
  int result;

  result = valence_one(mol);
  return(1);
}


int 
phase4(ums_type *mol)
{
  int count,i,j;
  double bond_length;
  int flag;
  
  for (count = 1; count <= Atoms; count ++)
  {
    switch(Redo(count))  
    {
    case 1:
      for (i = 0; i < Valence(count); i++)
      {
	j = Connection(count,i);
	bond_length = distance(Point(count),Point(j));
	if ((bond_length <= V2_C2_C_CUTOFF) && 
	    (Type(j)[0] == 'C'))
	  strcpy(Type(count),"C2");
	else
	  if ((bond_length <= V2_C2_N_CUTOFF) && 
	      (Type(j)[0] == 'N'))
	    strcpy(Type(count),"C2");
      }
      for (i = 0; i < Valence(count); i++)
      {
	j = Connection(count,i);
	bond_length = distance(Point(count),Point(j));
	if ((bond_length > V2_C3_C_CUTOFF) && 
	    (Type(j)[0] == 'C'))
	  strcpy(Type(count),"C3");
	else
	  if ((bond_length > V2_C3_N_CUTOFF) && 
	      (Type(j)[0] == 'N'))
	    strcpy(Type(count),"C3");
	  else
	    if ((bond_length > V2_C3_O_CUTOFF) && 
		(Type(j)[0] == 'O'))
	      strcpy(Type(count),"C3");
      }
      break;
    case 2:
      for (i = 0; i < Valence(count); i++)
      {
	j = Connection(count,i);
	bond_length = distance(Point(count),Point(j));
	if ((bond_length <= V2_N2_C_CUTOFF) && 
	    (Type(j)[0] == 'C'))
	  strcpy(Type(count),"Npl");
	else
	  if ((bond_length <= V2_N2_N_CUTOFF) && 
	      (Type(j)[0] == 'N'))
	    strcpy(Type(count),"Npl");
      }
      break;
    case 3:
      {
	flag = 0;
	for (i = 0; i < Valence(count); i++)
	{
	  j = Connection(count,i);
	  bond_length = distance(Point(count),Point(j));
	  if ((bond_length <= V2_C2_C_CUTOFF) && 
	      (Type(j)[0] == 'C'))
	  {
	    strcpy(Type(count),"C2");
	    flag = 1;
	  }
	  else
	    if ((bond_length <= V2_C2_N_CUTOFF) && 
		(Type(j)[0] == 'N'))
	    {
	      strcpy(Type(count),"C2");
	      flag = 1;
	    }
	}
	if (flag == 0)
	  for (i = 0; i < Valence(count); i++)
	  {
	    j = Connection(count,i);
	    bond_length = distance(Point(count),Point(j));
	    if ((bond_length > V2_C3_C_CUTOFF) && 
		(Type(j)[0] == 'C'))
	    {
	      strcpy(Type(count),"C3");
	      flag = 1;
	    }
	    else
	      if ((bond_length > V2_C3_N_CUTOFF) && 
		  (Type(j)[0] == 'N'))
	      {
		strcpy(Type(count),"C3");
		flag = 1;
	      }
	      else
		if ((bond_length > V2_C3_O_CUTOFF) && 
		    (Type(j)[0] == 'O'))
		{
		  strcpy(Type(count),"C3");
		  flag = 1;
		}
		else
		  if (flag == 0)
		    if ((bond_length > GEN_C3_C_CUTOFF) && 
		      (Type(j)[0] == 'C'))
		    {
		      strcpy(Type(count),"C3");
		      flag = 1;
		    }
	  }
      }
      break;
    }
  }
  return(1);
}

int 
phase5(ums_type *mol)
{
  int count,i, j;
  int flag;
  
  for (count = 1; count <= Atoms; count ++)
  {
    if (strcmp(Type(count),"C2") == 0)
    {
      flag = 0;
      for (i = 0; i < Valence(count); i++)
      {
	j = Connection(count,i);
	if ((strstr("C3   DC    HC   N3   N3+   O3   ",Type(j)) == NULL) &&
	    (strstr("Pac   Sac   Sox  C1   S3    Cac  ",Type(j)) == NULL))
	  flag = 1;
      }
      if (flag == 0) 
	strcpy(Type(count),"C3");
    }
  }
  return(1);
}


int 
phase6(ums_type *mol)
{
  int i,j,k,l,m,n;
  int conn;
  int no_plus;
  int protonated;

  for (i = 1; i <= Atoms; i ++)
  {
    no_plus = 1;
    protonated = TRUE;
    if (EQ(Type(i),"N3"))
    {
      for (j = 0; j < Valence(i); j++)
      {
	conn = Connection(i,j);

	/* If an unsaturated atom is attached to this nitrogen then it should be Npl */
	if ((Valence(i) == 2) && (IsUnsatType(Type(conn))))
	{
	  protonated = FALSE;
	  strcpy(Type(i),"Npl");
	  break;
	}

	/* If the attached atom is something other that C3, H or D the nitrogen is
	   not positively charged */
	if (NOTEQ(Type(conn),"C3") && (Atomic_number(conn) != 1))
	{
	  protonated = FALSE;
	}
      }
      if (protonated)
	strcpy(Type(i),"N3+");
    }
    /* look for guanadinium nitrogens */

    else 
      if (EQ(Type(i),"C2"))
      {
	/* First see if we have an sp2 carbon surrounded by 3 sp2
	   nitrogens */

	m = 0;
	for (j= 0; j < Valence(i); j++)
	{
	  k = Connection(i,j);
	  if ((EQ(Type(k),"Npl")) || (EQ(Type(k),"N2")) || (EQ(Type(k),"Ng+")))
	    m++;
	}
	if (m == 3)  
	{
	  strcpy(Type(i),"C+");
	  for (j= 0; j < Valence(i); j++)
	  {
	    k = Connection(i,j);
	    strcpy(Type(k),"Ng+");
	  }
	}
#if 0
	/* If we have 3 planar nitrogens connected to the sp3 carbon 
	   it is possible that this is a guanadinium group.  Now check
	   each of the nitrogens to see if they are connected to either
	   a C2 or Npl  */
	
	if (m == 3)  
	{
	  no_plus = FALSE;
	  for (j = 0; j < Valence(i); j++)
	  {
	    k = Connection(i,j);
	    if ((EQ(Type(k),"Npl")) || (EQ(Type(k),"N2")))
	    {
	      strcpy(Type(k),"Ng+");
	      for (l = 0; l < Valence(k); l++)
	      {
		n = Connection(k,l);
		if (n != i)
		{
		  if (EQ(Type(n),"C2") || EQ(Type(n),"Npl"))
		  {
		    no_plus = TRUE; 
		    break;
		  }
		}
	      }
	    }
	  }
	}
	if (no_plus == 1) 
	  for (j = 0; j < Valence(i); j++)
	  {
	    k = Connection(i,j);
	    if (strcmp(Type(k),"Ng+") == 0)
	      strcpy(Type(k),"Npl");
	  }
#endif
      }
      else 
	if (strcmp(Type(i),"Cac") == 0)
	  for (j = 0; j < Valence(i); j++)
	  {
	    k = Connection(i,j);
	    if ((strncmp(Type(k),"O",1) == 0) &&
		 (count_heavy_atoms(mol,k) == 1))
	      strcpy(Type(k),"O-");
	  }
  }
  return(1);
}


int 
type_hydrogens(ums_type *mol)
{
  int i, conn;
  
  for (i = 1; i <= Atoms; i++)
  {
    if (Type(i)[0] == 'H')
    {
      strcpy(Type(i),"H");
      conn = Connection(i,0);
      if (Type(conn)[0] == 'C')
      {
	strcpy(Type(i),"HC");
      }
    }
  }
  return(1);
}


int 
valence_four(ums_type *mol)
{
  int count;
  
  for (count = 1; count <= Atoms; count ++)
  {
    if ((Valence(count) == 4) && (IsOrganic(count)))
    {
      switch(Type(count)[0])
      {
      case 'C':
	if ((strcmp(Type(count),"C") == 0))
	  strcpy(Type(count),"C3");
	break;
      case 'N':
	if (count_free_ox(mol,count) >= 1)
	  strcpy(Type(count),"Nox");
	else
	  strcpy(Type(count),"N3+");
	break;
      case 'P':
	if (strlen(Type(count)) == 1)
	{
	  if (count_free_ox(mol,count) >= 2)
	    strcpy(Type(count),"Pac");
	  else
	    if (count_free_ox(mol,count) == 1)
	      strcpy(Type(count),"Pox");
	    else
	      strcpy(Type(count),"P3+");
	}
	break;
      case 'S':
	if (strcmp(Type(count),"S") == 0)
	{
	  if (count_free_ox(mol,count) >= 3)
	    strcpy(Type(count),"Sac");
	  else
	    if (count_free_ox(mol,count) >= 1)
	      strcpy(Type(count),"Sox");
	    else
	      strcpy(Type(count),"S");
	}
	break;
      case 'B':
	if (count_free_ox(mol,count) >= 3)
	  strcpy(Type(count),"Bac");
	if (count_free_ox(mol,count) >= 1)
	  strcpy(Type(count),"Box");
	else
	  strcpy(Type(count),"B");
	break;
      }
    }
  }
  return(1);
}

int 
valence_three(ums_type *mol)
{
  int count;
  int k,l,m;
  double angle1,angle2,angle3,avg_angle;
  
  for (count = 1; count <= Atoms; count ++)
  {  
    if ((Valence(count) == 3) && (IsOrganic(count)))
    {
      k = Connection(count,0);
      l = Connection(count,1);
      m = Connection(count,2);
      
      angle1 = bond_angle(Point(k),
			  Point(count),
			  Point(l));
      angle2 = bond_angle(Point(k),
			  Point(count),
			  Point(m));
      angle3 = bond_angle(Point(l),
			  Point(count),
			  Point(m));
      avg_angle = (angle1 + angle2 + angle3)/3;

      switch(Type(count)[0])
      {
      case 'C':
	if (avg_angle < SP3_MAX) 
	  strcpy(Type(count),"C3");
	else
	  if (count_free_ox(mol,count) >= 2)
	    strcpy(Type(count),"Cac");      
	  else 
	    strcpy(Type(count),"C2");      
	break;
      case 'N':
	if (avg_angle < SP3_MAX) 
	  strcpy(Type(count),"N3");
	else
	  if (count_free_ox(mol,count) >= 2)
	    strcpy(Type(count),"Ntr");      
	  else 
	    strcpy(Type(count),"Npl");      
	break;
      case 'B':
	if (count_free_ox(mol,count) >= 1)
	  strcpy(Type(count),"Box");      
	else
	  strcpy(Type(count),"B");      
	break;
      case 'S':
	if (strcmp(Type(count),"S") == 0)
	{
	  if (count_free_ox(mol,count) >= 1)
	    strcpy(Type(count),"Sox");      
	  else
	    strcpy(Type(count),"S3+");      
	}
	break;
      }
    }
  }
  return(1);
}


int 
valence_two(ums_type *mol)
{
  int count;
  int k,l;
  double angle1;
  
  for (count = 1; count <= Atoms; count ++)
  {
    if ((Valence(count) == 2) && (IsOrganic(count)))
    {
      k = Connection(count,0);
      l = Connection(count,1);  
      angle1 = bond_angle(Point(k),
			  Point(count),
			  Point(l));

      switch(Type(count)[0])
      {
      case 'C':
	if (strcmp(Type(count),"C") == 0)
	{
	  if (angle1 < SP3_MAX) 
	  {
	    strcpy(Type(count),"C3");
	    Redo(count) = 1;
	  }
	  else
	    if (angle1 < SP_MIN) 
	    {
	      strcpy(Type(count),"C2");
	      if (angle1 < MAY_BE_SP2)
		Redo(count) = 3;
	    }
	    else 
	      strcpy(Type(count),"C1");
	}
	break;
      case 'N':
	if (angle1 <= SP3_MAX) 
	{
	  strcpy(Type(count),"N3");
	  Redo(count) = 2;
	}
	else
	  if (angle1 <= SP_MIN) 
	  {
	    strcpy(Type(count),"Npl");
	  }
	  else 
	    strcpy(Type(count),"N1");
	break;
      case 'O':
	    strcpy(Type(count),"O3");
	break;
      case 'S':
	if (strcmp(Type(count),"S") == 0)
	  strcpy(Type(count),"S3");
	break;
      }
    }
  }
  return(1);
}      
  

int 
valence_one(ums_type *mol)
{
  int count;  
  int k;
  double bond_length;
  
  for (count = 1; count <= Atoms; count ++)
  {
    k = Connection(count,0);
    bond_length = distance(Point(count),Point(k));
    
    if ((Valence(count) == 1) && (IsOrganic(count)))
      switch(Type(count)[0])
      {
      case 'C':
	if (strcmp(Type(count),"C") == 0)
	{
	  if ((strncmp(Type(k),"C1",2) == 0) 
	    && (bond_length <= V1_C1_C1_CUTOFF))
	    strcpy(Type(count),"C1");
	  else
	    if ((strncmp(Type(k),"C",1) == 0)  
	      && (bond_length <= V1_C2_C_CUTOFF))
	      strcpy(Type(count),"C2");
	    else
	      strcpy(Type(count),"C3");
	}
	if (strncmp(Type(k),"N",1) == 0)
	{
	  if (bond_length <= V1_C2_N_CUTOFF)
	    strcpy(Type(count),"C2");
	  else
	    strcpy(Type(count),"C3");
	}
	break;
      case 'N':
	if (strcmp(Type(count),"N") == 0)
	  if ((strncmp(Type(k),"C1",2) == 0)
	  && (bond_length <= V1_N1_C1_CUTOFF))
	      strcpy(Type(count),"N1");
	  else
	    if (((strncmp(Type(k),"C2",2) == 0) ||
		(strncmp(Type(k),"C3",2) == 0))
	      && ((bond_length > V1_N3_C_CUTOFF)))
		strcpy(Type(count),"N3");
	    else
	      if (((strncmp(Type(k),"N3",2) == 0))
		&& ((bond_length > V1_N3_N3_CUTOFF)))
		  strcpy(Type(count),"N3");
	      else
		if (((strncmp(Type(k),"Npl",3) == 0))
		  && (bond_length > V1_N3_N2_CUTOFF))
		    strcpy(Type(count),"N3");
		else
		    strcpy(Type(count),"Npl");
	break;
      case 'O':
	if (strcmp(Type(count),"O") == 0)
	  if (strstr("Cac  Pac  Sac  Ntr  ",Type(k)) != NULL)
	    strcpy(Type(count),"O-");
	  else
	    if (strstr("Nox  Pox  Sox  ",Type(k)) != NULL)
	      strcpy(Type(count),"O2");
	    else
	      if ((Type(k)[0] == 'C')
		  && (bond_length <= V1_O2_C2_CUTOFF))
	      {
		strcpy(Type(count),"O2");
		strcpy(Type(k),"C2");
		Redo(k) = 0;
	      }
	      else
		if ((strcmp(Type(k),"As") == 0)
		    && (bond_length <= V1_O2_AS_CUTOFF))
		  strcpy(Type(count),"O2");
		else
		  strcpy(Type(count),"O3");
	break;
      case 'S':
	if (strcmp(Type(count),"S") == 0)
	if ((strncmp(Type(k),"P",1) == 0))
	  strcpy(Type(count),"S2");
	else
	  if ((strncmp(Type(k),"C",1) == 0)
	      && (bond_length <= V1_S2_C2_CUTOFF))
	  {
	    strcpy(Type(count),"S2");
	    strcpy(Type(k),"C2");
	    Redo(count) = 0;
	  }
	  else
	    if ((strcmp(Type(k),"As") == 0)
		&& (bond_length <= V1_S2_AS_CUTOFF))
	      strcpy(Type(count),"S2");
	    else
	      strcpy(Type(count),"S3");
	break;
      }
  }
  return(1);
}


double 
bond_angle(coord_type a, coord_type b, coord_type c)
{
  double angle;
  double dist;

  double cos_theta;

  dist = distance(a,b) * distance(b,c);
  cos_theta = ((a.x - b.x) * (c.x - b.x) + (a.y - b.y) * (c.y - b.y) +
	       (a.z - b.z) * (c.z - b.z))/dist;
  if (cos_theta  + 1.0 < 0.0001) 
    angle = 180.0;
  else
    angle = (acos(cos_theta)) * RAD_TO_DEG;
  return(angle);
}



	  
int 
count_heavy_atoms(ums_type *mol, int atom_number)
{
  int count;
  int H_count = 0;
  int bonded_atom;
    
  for (count = 0; count < Valence(atom_number); count ++)
  {
    bonded_atom = (Connection(atom_number,count));
    if (Type(bonded_atom)[0] == 'H')
    {
      H_count ++;
    }
  }
  return( Valence(atom_number) - H_count);
}


int 
count_free_ox(ums_type *mol, int atom_number)
{ 
  int count;
  int bonded_atom;
  int free_O_count = 0;

  for (count = 0; count < Valence(atom_number); count ++)
  {
    bonded_atom = (Connection(atom_number,count));
    if ((Type(bonded_atom)[0] == 'O') &&
	(count_heavy_atoms(mol,bonded_atom) == 1))
    {
      free_O_count ++;
    }
  }
  return(free_O_count);
}

void fix_carboxylates(ums_type *mol)
{
  int i;
  
  for (i = 1; i < Atoms; i++)
  {
    if EQ(Type(i),"O-")
    {
      switch(Valence(i))
      {
      case 1 :
	strcpy(Type(i),"O2");
	break;
      default :
	strcpy(Type(i),"O3");
	break;
      }
    }
  }
}

void check_for_amides(ums_type *mol)
{
  int i,j,conn;
  
  for (i = 1; i <= Atoms; i++)
  {
    if (EQ(Type(i),"Npl"))
      for (j = 0; j < Valence(i); j++)
      {
	conn = Connection(i,j);
	
	if (EQ(Type(conn),"Cac") || EQ(Type(conn),"Sox") || EQ(Type(conn),"So2"))
	{
	  strcpy(Type(i),"Nam");
	  break;
	}

	if (EQ(Type(conn),"C2"))
	  if (check_for_carbonyl(mol,Connection(i,j)) == 3)
	  {
	    strcpy(Type(i),"Nam");
	    break;
	  }
      }
  }
}










