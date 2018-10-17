/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : rdxyz.c
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : check the integrity of a structure
******/


#include "bbltyp.h"
    
int dummy_check(ums_type *mol)
{
  int i;
  
  for (i = 1; i <= Atoms; i++)
  {
    if (EQ(Type(i),"X"))
      return(FALSE);
    if ((X(i) > 10000.0) || (X(i) < -10000.0))
      return(FALSE);
    if ((Y(i) > 10000.0) || (Y(i) < -10000.0))
      return(FALSE);
    if ((Z(i) > 10000.0) || (Z(i) < -10000.0))
      return(FALSE);
  }
  return(TRUE);
}


