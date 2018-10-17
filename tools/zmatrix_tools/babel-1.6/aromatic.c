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

FILE : aromatic.c
AUTHOR(S) : Pat Walters
DATE : 2-96
PURPOSE : routines to find aromatic atom types and aromatic bonds
******/

#include "bbltyp.h"

#define AtomIsAromatic(x)  bit_is_on(info->arom_atms,x)

void find_aromatic_atoms(ums_type *mol)
{
  ring_struct rings;
  ring_info info;
  int i;

  find_SSSR(mol,&rings);
  setup_ring_info(mol,&rings,&info);
  for (i = 1; i <= Atoms; i++)
  {
    if (bit_is_on(info.arom_atms,i)) /* atom is aromatic */
    {
      switch(Atomic_number(i))
      {
      case 6 :
	strcpy(Type(i),"Car");
	break;
      case 7 :
	strcpy(Type(i),"Nar");
	break;
      }
    }
  }    

  for (i = 0; i < Bonds; i++)
  {
    if (bond_is_aromatic(Start(i),End(i),&info))
      Bond_order(i) = 5;
  }

  cleanup_ring_info(mol,&rings,&info);
  cleanup_rings(&rings);
}


/*------------------------------------------------
FUNCTION : count_heavy_atm_bonds
PURPOSE : count the number of bonds to non-hydrogen
atoms.  
-------------------------------------------------*/

int count_arom_atm_bonds(ums_type *mol, int atm, ring_info *info, int which_ring)
{
  int hvy_bonds = 0;
  int i, conn;
  
  for (i = 0; i < Valence(atm); i++)
  {
    conn = Connection(atm,i);
    if (in_same_ring(atm,conn,info,which_ring))
    { 
      if (Atomic_number(conn) != 1)
	hvy_bonds += BO(atm,i);  /* add bond order of bonds to non hydrogen atoms */
    }
  }
  return(hvy_bonds);  
}


void setup_ring_info(ums_type *mol, ring_struct *rings, ring_info *info)
{
  int i,j;

  info->num = NumRings;
  info->ring_atms = init_set_minbits(Atoms + 1);
  info->arom_atms = init_set_minbits(Atoms + 1);
  info->arom_rings = init_set_minbits(NumRings + 1);
  info->rings = (set_type **)malloc((Atoms + 1) * sizeof(set_type *));


  for (i = 1; i <= Atoms; i++)
  {
    info->rings[i] = init_set_minbits(NumRings + 1);
  }

  for (i = 0; i < NumRings; i++)
  {
    for (j = 0; j < RingSize(i); j++)
    {
      biton(info->rings[RingAtom(i,j)],i);
      biton(info->ring_atms,RingAtom(i,j));
    }
  }

  find_aromatic_rings(mol,rings,info);
  find_aromatic_rings2(mol,rings,info); 
}

void cleanup_ring_info(ums_type *mol, ring_struct *rings, ring_info *info)
{
  int i;
  
  free_set(info->ring_atms);
  free_set(info->arom_atms);
  free_set(info->arom_rings);
  
  for (i = 1; i <= Atoms; i++)
    free_set(info->rings[i]);
  
  free(info->rings);
}    

void print_ring_info(ums_type *mol, ring_info *info)
{
  char num_str[5];
  int i;

  printf("there are %d rings \n",info->num);
  setprint(info->arom_atms,"  arom_atms");
  setprint(info->arom_rings,"arom_rings");
  for (i = 1; i <= Atoms; i++)
  {
    sprintf(num_str,"%4d",i);
    setprint(info->rings[i],num_str);
  }
}

void find_aromatic_rings2(ums_type *mol, ring_struct *rings, ring_info *info)
{
  int i,j,is_aromatic,hyb;
  char hyb_str[5];
  
  for (i = 0; i < NumRings; i++)
  {
    is_aromatic = FALSE;
    if ((RingSize(i) == 5) || (RingSize(i) == 6))
    {
      is_aromatic = TRUE;
      for (j = 0; j < RingSize(i); j++)
      {
	/* check to see if the atom is sp2 
	   if not, the ring isn't aromatic */

	get_output_type(i,"HYB",Type(RingAtom(i,j)),hyb_str,all_caps);
	hyb = atoi(hyb_str);
	if (hyb != 2)
	{
	  is_aromatic = FALSE;
	  break;
	}
	
	/* check to see if the atom has 3 bonds to
	   heavy atoms in the same ring.
	   If not, the ring isn't aromatic */

	if (AtomIsAromatic(RingAtom(i,j)) == FALSE)
	{
	  if (count_arom_atm_bonds(mol,RingAtom(i,j),info, i) != 3)
	  {
	    is_aromatic = FALSE;
	    break;
	  }
	}

      }
    }
    
    if (is_aromatic)
    {
      biton(info->arom_rings,i);
      for (j = 0; j < RingSize(i); j++)
	biton(info->arom_atms,RingAtom(i,j));
    }
  }
}



void find_aromatic_rings(ums_type *mol, ring_struct *rings, ring_info *info)
{
  int i,j,is_aromatic,hyb;
  char hyb_str[5];
  
  for (i = 0; i < NumRings; i++)
  {
    is_aromatic = FALSE;
    if ((RingSize(i) == 5) || (RingSize(i) == 6))
    {
      is_aromatic = TRUE;
      for (j = 0; j < RingSize(i); j++)
      {
	/* check to see if the atom is sp2 
	   if not, the ring isn't aromatic */

	get_output_type(i,"HYB",Type(RingAtom(i,j)),hyb_str,all_caps);
	hyb = atoi(hyb_str);
	if (hyb != 2)
	{
	  is_aromatic = FALSE;
	  break;
	}

	/* check to see if the atom has 3 bonds to
	   heavy atoms in the same ring.
	   If not, the ring isn't aromatic */

	if (count_arom_atm_bonds(mol,RingAtom(i,j),info, i) != 3)
	{
	  is_aromatic = FALSE;
	  break;
	}
      }
    }
    
    if (is_aromatic)
    {
      biton(info->arom_rings,i);
      for (j = 0; j < RingSize(i); j++)
	biton(info->arom_atms,RingAtom(i,j));
    }
  }
}




/*--------------------------------------------------------
FUNCTION - bond_is_aromatic

Check to see if a bond is aromatic.  Bond is considered aromatic
if both atoms are aromatic and both atoms are in the same
ring.
--------------------------------------------------------*/

int bond_is_aromatic(int a, int b, ring_info *info)
{
  set_type *common;
  int next, result = FALSE;
  
  /* Check to see if both atoms in the bond are aromatic */
  
  if ((!AtomIsAromatic(a)) || (!AtomIsAromatic(b)))
    return(FALSE);
  
  /* If both atoms are in the same ring the bond is aromatic */
  
  if (in_same_ring(a,b,info,-1))
    return(TRUE);

  return(FALSE);
}


int in_same_ring(int a, int b, ring_info *info, int which_ring)
{
  int result = FALSE;
  
  set_type *common, *aromatic_common;

  common = init_set_minbits(info->num);
  aromatic_common = init_set_minbits(info->num);
  setand(info->rings[a],info->rings[b],common);  

  if ((which_ring >= 0) && (bit_is_on(common,which_ring)))
    result = TRUE;
  else
    if ((which_ring < 0) && (setcount(common) > 0))
    {
      setand(common,info->arom_rings,aromatic_common);
      if (setcount(aromatic_common) > 0)
	result = TRUE;
    }
  
  
  free_set(common);
  free_set(aromatic_common);      
  return(result);
}


int count_ortho_substituents(ums_type *mol, int start, int end)
{
  int i, conn, ortho_count = 2;
  
  for (i = 0; i < Valence(start); i++)
  {
    conn = Connection(start,i);
    if (conn != end)
    {
      ortho_count -= count_bonded_hydrogens(mol,conn);
    }
  }

  return(ortho_count);
}

int count_bonded_hydrogens(ums_type *mol, int atm)
{
  int i, h_count = 0, conn;
  
  for (i = 0; i < Valence(atm); i++)
  {
    conn = Connection(atm,i);
    if (Atomic_number(conn) == 1)
      h_count++;
  }
  return(h_count);
}

