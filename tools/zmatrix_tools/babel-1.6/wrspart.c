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

FILE : wrspart.c
AUTHOR(S) : Pat Walters
DATE : 2-94
PURPOSE : routines to write a Spartan file
******/

#include "bbltyp.h"


int 
write_spartan(FILE *file1, ums_type *mol)
{ 
  int i;
  char ele[5],temp1[5];
  int atnum;
  
  if (NOTEQ(OutputKeywords,NOKEY))
    fprintf(file1,"%s\n",OutputKeywords);
  else
    fprintf(file1,"\n");
  fprintf(file1,"%s\n0 1\n",Title);
  for(i = 1;i <= Atoms; i++)
  {
    get_element_type(mol,i,ele);
    atnum = get_atomic_number(ele);
    fprintf(file1,"%2d%14.9f%14.9f%14.9f\n",
	    atnum,
	    X(i),
	    Y(i),
	    Z(i));
  }
  fprintf(file1,"ENDCART\nPAIRING\nENDPAIR\nHESSIAN\n");
  for (i = 1; i <= Atoms; i++)
  {
    get_output_type(i,"MOL",Type(i),temp1,dummy);
    fprintf(file1,"%5d",atoi(temp1));
    if ((i % 12) == 0)
      fprintf(file1,"\n");
  }
  if ((Atoms % 12) != 0)
    fprintf(file1,"\n");
  for(i = 0; i < Bonds; i++)
  {
    fprintf(file1,"%5d%5d%5d\n",
	    Start(i),
	    End(i),
	    Bond_order(i));
  }
  fprintf(file1,"ENDHESS\n");
  return(TRUE);
}


  
      









