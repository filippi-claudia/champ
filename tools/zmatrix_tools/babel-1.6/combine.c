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

FILE : combine.c
AUTHOR(S) : Pat Walters
DATE : 11-92
PURPOSE : routines to combine the contents of 2 UMS structures

******/



#include "bbltyp.h"
static warning wstr;

int 
borrow_types(ums_type *mol1, ums_type *mol2)
{
  int count1, count2;
  int replaced;
  
  for (count1 = 1; count1 <= mol1->num_atoms; count1 ++)
  {
    replaced = FALSE;
    for (count2 = 1; count2 <= mol2->num_atoms; count2 ++)
    {
      if ((mol1->atoms[count1].point.x == mol2->atoms[count2].point.x) &&
	  (mol1->atoms[count1].point.y == mol2->atoms[count2].point.y) &&
	  (mol1->atoms[count1].point.z == mol2->atoms[count2].point.z))
      {
	strcpy(mol1->atoms[count1].type,mol2->atoms[count2].type);
	replaced = TRUE;
	break;
      }
    }
    if (replaced == FALSE)
      sprintf(wstr,"Could not replace type for atom %d",count1);
    show_warning(wstr);
  }
  return(TRUE);
}

      
