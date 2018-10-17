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

FILE : delatms.c
AUTHOR(S) : Matt Stahl
DATE : 9-8-94
PURPOSE : remove Oxygens with valence 0 from a UMS

******/

#include "bbltyp.h"

ums_type *delete_water(ums_type *mol)
 {  
  int heavy_count;

  heavy_count = tag_waters(mol);
  return(build_new_ums(mol,heavy_count));
}

int tag_waters(ums_type *mol)
{
  int i;
  int j = 0;

  for (i = 1;i <= Atoms;i ++)
  {   
    if (is_water(mol,i))
    {
      Redo(i) = 0;
    }
    else
    {
      j ++;
      Redo(i) = j;
    }
  }
  return(j);
}

int is_water(ums_type *mol, int i)
{
  if (EQ(ResName(i),"HOH"))
    return(TRUE);
  if (EQ(ResName(i),"H2O"))
    return(TRUE);
  if (EQ(ResName(i),"WAT"))
    return(TRUE);
  if (EQ(ResName(i),"OH2"))
    return(TRUE);
  return(FALSE);
}



