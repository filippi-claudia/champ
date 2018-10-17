/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : ringutil.c
AUTHOR(S) : Matt Stahl, Pat Walters
DATE : 2-10-93
PURPOSE : utilites for ring perception
******/


#include "bbltyp.h"

int find_SSSR(ums_type *mol, ring_struct *rings)
{
  int i;
  ums_type *new_u,*head,*tmp;
  ring_struct tmp_rings;
  int count = 0;
  
  rings->ring_list = NULL;
  rings->count = 0;

  tmp = (ums_type *)malloc(sizeof(ums_type));
  if (!tmp)
    fatal_error("Unable to allocate memory for temporary UMS");
  
  head = new_u = dissect_ums(mol);

  while (new_u)
  {
    if (new_u->num_atoms > 2)
    {
      tmp->num_atoms = new_u->num_atoms;
      tmp->num_bonds = new_u->num_atoms;
      
      initialize_ums(&tmp);
      
      for (i = 1;i <= new_u->num_atoms;i++)
	tmp->atoms[i].redo = new_u->atoms[i].redo;
      
      tmp_rings.count = 0;

      find_rings(new_u,&tmp_rings);

      add_rings_to_list(rings,&tmp_rings,tmp);
      
      if (tmp_rings.ring_list)
	free(tmp_rings.ring_list);
      
      release_ums(tmp);
    }
    count++;
    new_u = new_u->next;
  }

  if (tmp)
    free(tmp);

  new_u = head;
  while (new_u)
  {
    tmp = new_u->next;
    release_ums(new_u);
    free(new_u);
    new_u = tmp;
  }

  if (new_u)
    free(new_u);
  return(TRUE);
}
ums_type *dissect_ums(ums_type *mol)
{
  register int i,j;
  ums_type *head,*ptr;
  set_type *curr = init_set_minbits(Atoms);
  set_type *next = init_set_minbits(Atoms);
  set_type *used = init_set_minbits(Atoms);
  set_type *mask = init_set_minbits(Atoms);

  head = NULL;
  while (setcount(used) < Atoms)
  {
    setclear(curr);
    for (i = 1;i <= Atoms;i++)
       if (!bit_is_on(used,i))
       {
	  biton(curr,i);
	  break;
       }

    setclear(mask);
    setor(curr,mask,mask);
    setor(curr,used,used);

    do
    {
      setclear(next);
      for (i = NextBit(curr,0);i != -1;i = NextBit(curr,i))
	 for (j = 0;j < Valence(i);j++)
	    if (!bit_is_on(used,Connection(i,j)))
	       biton(next,Connection(i,j));
       
      setcopy(curr,next);
      setor(curr,used,used);
      setor(curr,mask,mask);
    } while (setcount(curr) > 0);

    ptr = set_to_ums(mol,mask);
    head = add_ums_to_list(ptr,head);
  }

  free_set(curr);
  free_set(next);
  free_set(used);
  free_set(mask);

  return(head);
}

void add_rings_to_list(ring_struct *rings,ring_struct *tmp,ums_type *mol)
{
  register int i,j,count,bit;
  int last;

  if (rings->ring_list == NULL)
  {
    rings->ring_list = (path *)malloc(sizeof(path) * tmp->count);
    if ((tmp->count > 0) && (rings->ring_list == NULL))
      fatal_error("Unable to allocate memory for ring structures\n");
    rings->count = 0;
  }
  else
  {
    rings->ring_list = 
      (path *)realloc(rings->ring_list,(sizeof(path) * 
					(rings->count + tmp->count)));
    if ((tmp->count > 0) && (rings->ring_list == NULL))
      fatal_error("Unable to allocate memory for ring structures");
  }
  count = rings->count;
  preserve_rings(rings,tmp,tmp->count);

  for (i = count;i < rings->count;i++) 
  {
    for (j = 0;j < rings->ring_list[i].length;j++)
    {
      rings->ring_list[i].path_atoms[j] = Redo(rings->ring_list[i].path_atoms[j]);
    }
      
    free_set(rings->ring_list[i].path_set);
    
    last = find_last_atom(&rings->ring_list[i]);
    rings->ring_list[i].path_set = init_set_minbits(last);

    for (j = 0;j < rings->ring_list[i].length;j++)
    {
    	bit = rings->ring_list[i].path_atoms[j];
      biton(rings->ring_list[i].path_set,bit);
    }

  }
}


ums_type *set_to_ums(ums_type *mol,set_type *set)
{
  int i,atom = 1;
  ums_type *new_u;
  
  new_u = (ums_type *)malloc(sizeof(ums_type));
  
  if (new_u == NULL)
  {
    fprintf(stderr,"Unable to allocate memory for temporary ums\n");
    exit(0);
  }
  
  new_u->num_atoms = setcount(set);
  new_u->num_bonds = 0;
  initialize_ums(&new_u);

  strcpy(new_u->title,Title);
  
  for (i = 1;i <= Atoms;i++)
    if (bit_is_on(set,i))
    {
      Redo(i) = atom;
      atom++;
    }
    else
      Redo(i) = 0;
  
  if (Atoms > 0)
  {
    for (i = 0;i < Bonds;i++)
    {
      if (bit_is_on(set,Start(i)) && bit_is_on(set,End(i)))
      {
	new_u->connections[new_u->num_bonds].start = Redo(Start(i));
	new_u->connections[new_u->num_bonds].end = Redo(End(i));
	new_u->connections[new_u->num_bonds].bond_order = Bond_order(i);
	new_u->atoms[Redo(Start(i))].redo = Start(i);
	new_u->atoms[Redo(End(i))].redo = End(i);
	new_u->num_bonds++;
      }
    }
    
    for (i = 1;i <= Atoms;i++)
    {
      if (Redo(i) != 0)
      {
	if (!isalpha(Type(i)[0]))
	  fatal_error("Something is wrong with the dissected ums");
	
	strcpy(new_u->atoms[Redo(i)].type,Type(i));
	new_u->atoms[Redo(i)].point.x = X(i);
	new_u->atoms[Redo(i)].point.y = Y(i);
	new_u->atoms[Redo(i)].point.z = Z(i);
	new_u->atoms[Redo(i)].charge = Charge(i);
	new_u->atoms[Redo(i)].dble = Double(i);
	new_u->atoms[Redo(i)].radius = Radius(i);
	new_u->atoms[Redo(i)].atomic_number = Atomic_number(i);
	new_u->atoms[Redo(i)].redo = i;
      }
    }
    dissect_connection_table(new_u);
  }
  return(new_u);
}

int find_last_atom(path *the_path)
{
  int i;
  int last = -1;
  
  for (i = 0; i < the_path->length; i++)
  {
    if (the_path->path_atoms[i] > last)
      last = the_path->path_atoms[i];
  }
  return(last);
}   

void cleanup_rings(ring_struct *rings)
{
  int i;
  
  for (i = 0; i < rings->count; i++)
  {
    free_set(rings->ring_list[i].path_set);
    if (rings->ring_list[i].path_atoms)
      free(rings->ring_list[i].path_atoms);
  }
  if (rings->ring_list)
    free(rings->ring_list);
}


void preserve_rings(ring_struct *good,ring_struct *bad,int count) 
{
  int i,j,k,*tmp;
  set_type *set;
  int duplicate;
  
  j = good->count;
  for (i = 0; i < count; i++)
  {
    duplicate = FALSE;
    for (k = 0;k < j;k++)
      if (setcmp(good->ring_list[k].path_set,bad->ring_list[i].path_set))
	duplicate = TRUE;

    if (!bad->ring_list[i].bogus)
    {
      tmp = bad->ring_list[i].path_atoms;
      set = bad->ring_list[i].path_set;
      
      bad->ring_list[i].path_atoms = NULL;
      bad->ring_list[i].path_set = NULL;
      
      good->ring_list[j] = bad->ring_list[i];
      good->ring_list[j].path_atoms = tmp;
      good->ring_list[j].path_set = set;
      j++;
    }
    else
    {
      if (bad->ring_list[i].path_atoms)
	free(bad->ring_list[i].path_atoms);
      free_set(bad->ring_list[i].path_set);
    }
  }      
  
  good->count = j;
}




