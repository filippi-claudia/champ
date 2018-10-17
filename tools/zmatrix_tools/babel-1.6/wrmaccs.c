/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
------------------------------------------------------------------------------

FILE : wrmaccs.c
AUTHOR(S) : Matt Stahl
DATE : 11-92
PURPOSE : Routines to write the maccs format as defined by mdl
******/

#include "bbltyp.h"

int 
write_maccs(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  int result;

  fprintf(file1,"%3d%3d\n",Atoms,Bonds);

  for(i = 1;i <= Atoms; i++)
  {

    switch ((int) Charge(i))
    {
    case 3:
      Charge(i) = 1.0;
      break;
    case 1:
      Charge(i) = 3.0;
      break;
    case -1:
      Charge(i) = 5.0;
      break;
    case -2:
      Charge(i) = 6.0;
      break;
    case -3:
      Charge(i) = 7.0;
      break;
    }

    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1,"%10.4f%10.4f%10.4f %-3s 0%3d  0  0  0\n",
	    X(i),
	    Y(i),
	    Z(i),
	    type_name,
	    (int)Charge(i));
  }

  for(i = 0;i < Bonds; i++)
    fprintf(file1,"%3d%3d%3d  0  0  0\n",Start(i),End(i),Bond_order(i));
  fprintf(file1,"$$$$\n");
  return(TRUE);
}


  
      









