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
FILE : rdhin.c
AUTHOR(S) : Pat Walters
DATE : 1-94
PURPOSE : routines to read a Hyperchem .HIN file
******/
#include "bbltyp.h"
#define DELIMS "\t\n "

int 
read_hin(FILE *file1, ums_type *mol)
{
  char the_line[300];
  int tokens;
  int i,j,k;
  int max;
  char temp[20];
  long pos = 0;
  int counter = 0;

  while (strstr(the_line,"mol") == NULL)
    fgets(the_line,sizeof(the_line),file1);
  Atoms = 0;
  if (strstr(the_line,"mol ") != NULL)
  {
    counter++;
    pos = ftell(file1);
    Atoms = 0;
    while (strstr(the_line,"endmol") == NULL)
    {
      fgets(the_line,sizeof(the_line),file1);
      if (strstr(the_line,"atom"))
      {
	Atoms++;
      }
    }
  }

  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  fseek(file1,pos,0);
  i = 1;
  strcpy(the_line,"");
  while (strstr(the_line,"endmol") == NULL)
  {
    fgets(the_line,sizeof(the_line),file1);
    if (strstr(the_line,"atom"))
    {
      UpdateProgress();
      tokens = count_tokens(the_line,DELIMS);
      get_token(Type(i),the_line,DELIMS,4);
      clean_atom_type(Type(i));
      get_token(temp,the_line,DELIMS,8);
      X(i) = atof(temp);
      get_token(temp,the_line,DELIMS,9);
      Y(i) = atof(temp);
      get_token(temp,the_line,DELIMS,10);
      Z(i) = atof(temp);
      get_token(temp,the_line,DELIMS,11);
      Valence(i) = (int) atof(temp);
      k = 0;
      max = 12 + 2 * Valence(i);
      for (j = 12; j < max; j+=2)
      {
	get_token(temp,the_line,DELIMS,j);
	Connection(i,k) = atoi(temp);
	get_token(temp,the_line,DELIMS,j+1);
	BO(i,k) = translate_hin_bond_order(temp[0]);
	k++;
      }
      i++;
    }
  }
  
  assign_types(mol);
  build_connection_table(mol);
  
  return(TRUE);
}


int translate_hin_bond_order(char c)
{
  int bo;
  
  switch(c)
  {
  case 's':
    bo = 1;
    break;
  case 'd':
    bo = 2;
    break;
  case 't':
    bo = 3;
    break;
  case 'a':
    bo = 5;
    break;
  default :
    printf("Can't translate .hin bond order %c\n",c);
  }
  return(bo);
}

    






   
    
    
    
	  


