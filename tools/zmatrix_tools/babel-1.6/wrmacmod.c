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

FILE : wrmacmod.c
AUTHOR(S) : Pat Walters
DATE : 2-93
PURPOSE : Routines to write a Macromodel file
******/


#include "bbltyp.h"

int 
  write_macromodel(FILE *file1, ums_type *mol)
{ 
  int i,j;
  int type_name;
  char temp_type[5];
  int result;

  fprintf(file1," %5d %6s      E = %7.3f KJ/mol\n",Atoms,Title,4.184*Energy);
  
  for(i = 1;i <= Atoms; i++)
  {
    if (Atomic_number(i) == 1)
      switch(Atomic_number(Connection(i,0)))
      {
      case '7' :
	type_name = 43;
	break;
      case '8' :
	type_name = 42;
	break;
      default :
	type_name = 41;
      }
    else
    {    
      result = get_output_type(i,"MMD",Type(i),temp_type,dummy);    
      type_name = atoi(temp_type);
    }
    fprintf(file1,"%4d",type_name);
    for (j = 0; j < 6; j ++)
    {
      fprintf(file1," %5d %1d",
	      Connection(i,j),BO(i,j));
    }
    fprintf(file1," %11.6f %11.6f %11.6f %5d %5d %8.5f %8.5f\n",
	    X(i),
	    Y(i),
	    Z(i),
	    0,
	    0,
	    Charge(i),
	    Charge(i));
  }
  return(TRUE);
}

int 
get_bond_order(ums_type *mol, int start,int end)
{
  int i;
  
  if ((start == 0) || (end == 0))
    return(0);
  for (i = 0; i < Bonds; i++)
    if (((start == Start(i)) && (end == End(i))) ||
	((start == End(i)) && (end == Start(i))))
      return(Bond_order(i));
  return(0);
}
  
  
  
      









