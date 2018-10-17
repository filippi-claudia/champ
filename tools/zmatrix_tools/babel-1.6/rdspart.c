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
FILE : rdspart.c
AUTHOR(S) : Pat Walters
DATE : 
PURPOSE : routines to read a Spartan Input file
******/

#include "bbltyp.h"

extern element_type *elements;

int 
read_spartan(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i = 1,j,k,l;
  int skip = -1;
  int atnum;
  char temp1[5],temp2[5];
  int start = 1, end;
  int column;

  Atoms = 0;
  while ((fgets(the_line,sizeof(the_line), file1) != NULL) &&
	 (strstr(the_line,"ENDCART") == NULL))
  {
    i++;
    if ((i > 4) && (count_tokens(the_line,"\n\t ") == 4))
    {
      if (skip < 0)
	skip = i - 1;
      Atoms++;
    }
  }
  printf("there are %d atoms\n",Atoms);
  initialize_ums(&mol);

  rewind(file1);
  for (i = 1; i < skip; i++)
    fgets(the_line,sizeof(the_line),file1);
  for (i = 1; i <= Atoms; i++)
  {
    fgets(the_line,sizeof(the_line),file1);
    sscanf(the_line,"%d %lf %lf %lf",&atnum,&X(i),&Y(i),&Z(i));
    strcpy(Type(i),elements[atnum].name);
  }

  while ((fgets(the_line,sizeof(the_line), file1) != NULL) &&
	 (strstr(the_line,"HESSIAN") == NULL));

  column = locate_input_type("MOL");

  l = 0;
  for (i = 0; i < ceil(Atoms/12.0); i++)
  {
    end = start + 12;
    if (end > Atoms)
      end = Atoms;
    fgets(the_line,sizeof(the_line),file1);
    k = 0;
    for (j = start; j <= end; j++)
    {
      l++;
      strcpy(temp2,Type(j));
      sscanf(&the_line[k],"%s",temp1);
      if (atoi(temp1) > 0)
	Atomic_number(j) = get_input_type(j,column,temp1,Type(j),dummy);    
      if (EQ(Type(i),"Du"))
	strcpy(Type(j),temp2);
      k+=5;
    }
    start+=12;
  }
  i = 0;
  while ((fgets(the_line,sizeof(the_line), file1) != NULL) &&
	 (strstr(the_line,"ENDHESS") == NULL))
    {
      sscanf(the_line,"%d %d %d",&Start(i),&End(i),&Bond_order(i));
      i++;
    }
  Bonds = i;
  dissect_connection_table(mol);
  assign_bond_order(mol);
  read_to_eof(file1);
  return(TRUE);
}







   
    
    
    
	  


