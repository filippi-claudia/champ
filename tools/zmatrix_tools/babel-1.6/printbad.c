/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
----------------------------------------------------------------------------

FILE : printbad.c
AUTHOR(S) : Pat Walters
DATE : 9-2-94
PURPOSE : print out bond distances for atoms which execeed max valence
******/

#include "bbltyp.h"

void print_bad_connections(ums_type *mol, int atom)
{
  int i;
  double dist;
  int conn;
  
  for (i = 0; i < Valence(atom); i++)
  {
    conn = Connection(atom,i);
    dist = distance(Point(atom),Point(conn));
    if (mol->residues != NULL)
      fprintf(stderr,"%4d - %4d %4s in residue %3s%3d %4s in residue %3s%3d distance = %9.3f\n",
	      atom,conn,AtmId(atom),ResName(atom),ResNum(atom),
	      AtmId(conn),ResName(conn),ResNum(conn),dist);
    else
      fprintf(stderr,"%4d - %4d %4s %4s distance = %9.3f\n",
	      atom,conn,Type(atom),Type(conn),dist);
  }
}

