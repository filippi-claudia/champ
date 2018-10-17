/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
-----------------------------------------------------------------------------
*/
#include "bbltyp.h"


int do_spline(ums_type *mol)
{		
  FILE *file1, *outfile;
  int end = FALSE, result = 0;
  int i;
  int in_count = 0;
  int out_count = 1;
  ums_type *mol2,*tmp;
  vect_type **vect;
  
  mol2 = (ums_type *)malloc(sizeof(ums_type));
  if (!mol2)
    fatal_error("Unable to allocate memory for ums");
  mol2->control = mol->control;
  
  file1 = open_read(InfileName);
  
  vect = NULL;
  while (!end)
  {
    if (in_count == 0)
    {
      result = ReaderFunction(file1,mol);
      result = ReaderFunction(file1,mol2);
      in_count = 2;
      vect = (vect_type **)malloc(sizeof(vect_type *) * (Atoms+1));
      for (i = 1;i <= Atoms;i++)
	vect[i] = (vect_type *)malloc(sizeof(vect_type));
    }
    else
    {
      tmp = mol;
      mol = mol2;
      mol2 = tmp;
      release_ums(mol2);
      result = ReaderFunction(file1,mol2);
      in_count++;
    }
    
    get_vectors(mol,mol2,vect,Increment);

    for (i = 0;i < Increment;i++)
    {
      generate_outfile_name(mol,out_count);
      outfile = open_write(OutfileName);
      do_outputs(outfile,mol);
      fclose(outfile);
      out_count++;
      add_step(mol,vect);
    }

    if (check_for_eof(file1))
      end = TRUE;
  }

  return(FALSE);
}

void
  get_vectors(ums_type *start,ums_type *end,vect_type *vect[],int increment)
{
  int i;
  
  if (increment == 0)
    increment = 1;
    
  for (i = 1;i <= start->num_atoms;i++)
  {
    (*vect[i]).x = start->atoms[i].point.x - end->atoms[i].point.x;
    (*vect[i]).y = start->atoms[i].point.y - end->atoms[i].point.y;
    (*vect[i]).z = start->atoms[i].point.z - end->atoms[i].point.z;
  }

  for (i = 1;i <= start->num_atoms;i++)
  {
    (*vect[i]).x  = (*vect[i]).x/(double) increment;
    (*vect[i]).y  = (*vect[i]).y/(double) increment;
    (*vect[i]).z  = (*vect[i]).z/(double) increment;
  }
}

void
  add_step(ums_type *mol,vect_type *vect[])
{
  int i;
  
  for (i = 1;i <= Atoms;i++)
  {
    X(i) += -(*vect[i]).x;
    Y(i) += -(*vect[i]).y;
    Z(i) += -(*vect[i]).z;
  }
}
