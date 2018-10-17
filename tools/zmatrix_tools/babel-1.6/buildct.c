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

FILE : buildct.c
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : Build a connection table from a UMS

******/



#include "bbltyp.h"

int 
build_connection_table(ums_type *mol)
{
  int i,j;
  int num_bonds = 0;

  for (i = MIN_ATOM; i < Atoms; i++)
  {
    for (j = 0; j < Valence(i); j++)
    {
      if (member(i,Connection(i,j),mol->connections,num_bonds) == 1)
      {
	Start(num_bonds) = i;
	End(num_bonds) = Connection(i,j);
	Bond_order(num_bonds) = BO(i,j);
	num_bonds++;
      }
    }
  }
  Bonds = num_bonds;
  return(TRUE);
}


int 
member(int first, int second, connect_type *connections, int length)
{
  int count;

  if (length == 0) return(1);
  for (count = 0; count < length; count++)
  {
    if ((first == connections[count].start) && (second == connections[count].end)) return(0);
    if ((second == connections[count].start) && (first == connections[count].end)) return(0);
  }
  return(1);
}



























