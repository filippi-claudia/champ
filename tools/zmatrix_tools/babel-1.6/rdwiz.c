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
FILE : rdwiz.c
AUTHOR(S) : Matthew Stahl
DATE : 5-94
PURPOSE : File format for use with Wizard
******/

#include "bbltyp.h"


int 
read_wizard(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE],title[BUFF_SIZE];
  int i;
  int result;
  int tokens;
  double strain;
  int column;
  char temp_type[10];
  
  fgets(the_line,sizeof(the_line), file1);
  strcpy(title,the_line);
  
  if (count_tokens(title," \n\t") == 0)
    strcpy(title,InfileName);
  strip_return(title);

  fgets(the_line,sizeof(the_line), file1);
  sscanf(the_line,"%lf",&strain);

  fgets(the_line,sizeof(the_line), file1);
  sscanf(the_line,"%d %d",&Atoms,&Bonds);

  ShowProgress(Atoms,"Reading Atoms");

  result = initialize_ums(&mol);
  Energy = strain;
  strcpy(Title,title);
  
  column = locate_input_type("INT");  
  for (i = MIN_ATOM; i <= Atoms;i ++)
    {
      UpdateProgress();
      fgets(the_line,sizeof(the_line), file1);
      tokens = count_tokens(the_line,"\t\n ");
      Valence(i) = tokens - 4;
      sscanf(the_line,"%s %lf %lf %lf",
	     temp_type,
	     &X(i),
	     &Y(i),
	     &Z(i));
      Atomic_number(i) = get_input_type(i,column,temp_type,Type(i),all_caps);    
    }

  for (i = 0;i < Bonds;i++)
  {
    fgets(the_line,sizeof(the_line), file1);
    sscanf(the_line,"%d %d %d",&Start(i),&End(i),&Bond_order(i));
  }

  dissect_connection_table(mol);
  OrderConnections(mol);

  return(TRUE);
}

void
OrderConnections(ums_type *mol)
{
  int i,j,k,temp,tempBO;

  for (i = 1;i <= Atoms;i++)
  {
    for (j = 0;j < Valence(i);j++)
    {
      if (EQn("H",Type(Connection(i,j)),1))
      {
	temp = Connection(i,j);
	tempBO = BO(i,j);
	
	for (k = j;k < (Valence(i) - 1);k++)
	{
	  Connection(i,k) = Connection(i,k+1);
	  BO(i,k) = BO(i,k+1);
	}
	Connection(i,k) = temp;
	BO(i,k) = tempBO;
      }
    }
  }
}






   
    
    
    
	  


