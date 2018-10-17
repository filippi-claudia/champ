/*-----------------------------------------------------------
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
-------------------------------------------------------------

FILE : umslist.c
AUTHOR(S) : Pat Walters
DATE : 1-96
PURPOSE : utilites for dealing with a list of ums structures
-------------------------------------------------------------*/

#include "bbltyp.h"

ums_type *add_ums_to_list(ums_type *new_node, ums_type *list)
{
  if (list == NULL)
  {
    list = new_node;
    return(list);
  }
  else 
  {
    if (list->num_atoms <= new_node->num_atoms)
    {
      new_node->next = list;
      list = new_node;
      return(list);
    }
    else
    {
      list->next = add_ums_to_list(new_node, list->next);
      return(list);
    }
  }
}

int count_ums_list(ums_type *list)
{
  int i = 1;
  
  if (list && !list->next)
    return(1);

  while (list->next)
  {
    i++;
    list = list->next;
  }
  return(i);
}

void cleanup_ums_lst(ums_type *base)
{
  ums_type *listp, *nextp;
  
  for (listp = base; listp != NULL; listp = nextp)
  {
    nextp = listp->next;
    release_ums(listp);
    if (listp)
      free(listp);
  }
}


ums_type *ums_list_from_file(FILE *fp,int (*reader)())
{
  ums_type *mol = NULL, *list = NULL;
  
  while (!check_for_eof(fp))
  {
    mol = (ums_type *)malloc(sizeof(ums_type));
    mol->next = NULL;
    reader(fp,mol);
    list = add_ums_to_list(mol,list);
  }
  return(list);
}

void ums_list_to_file(FILE *fp,ums_type *list,int (*writer)())
{
  while (list)
  {
    writer(fp,list);
    list = list->next;
  }
}

void show_ums_list(ums_type *list)
{
  while (list)
  {
/*
    write_xyz(stdout,list);
*/
    printf("%15s %10.3f\n",list->title,list->energy);
    list = list->next;
  }
}
