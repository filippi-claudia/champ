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

FILE : rdmm2in.c
AUTHOR(S) : Pat Walters
DATE : 1-93
PURPOSE : routines to read an mm2 input file, this is also the format of the
MM2 87 output file
	
******/

#include "bbltyp.h"

int 
read_mm2_input(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  char temp[25];
  int i,j;
  int METHOD,IBETA,LRPI,NCONN,NATTACH,LABEL,NDC,ICONT;
  int conn,attach;
  int tokens;
  int a,b;
  int column;

#ifdef AICHEM
  fgets(the_line,sizeof(the_line),file1);
#endif
  fgets(the_line,sizeof(the_line),file1);
  my_strncpy(temp,&the_line[60],1);
  METHOD = atoi(temp);
  my_strncpy(temp,&the_line[61],4);
  Atoms = atoi(temp);
  ShowProgress(Atoms,"Reading Atoms");

  initialize_ums(&mol);
  if (METHOD == 1)
  {
    fgets(the_line,sizeof(the_line),file1);
    my_strncpy(temp,&the_line[67],1);      
    IBETA = atoi(temp);
    my_strncpy(temp,&the_line[74],1);      
    LRPI = atoi(temp);
    my_strncpy(temp,&the_line[79],1);      
    ICONT = atoi(temp);
    if (ICONT == 1) /* CARD 1b */
      fgets(the_line,sizeof(the_line),file1);
    if (IBETA == 1) /* CARD 1c */
      fgets(the_line,sizeof(the_line),file1);
    if (LRPI != 0)  /* CARD 1d */
      fgets(the_line,sizeof(the_line),file1);
  }
  fgets(the_line,sizeof(the_line),file1);
  my_strncpy(temp,the_line,5);      
  NCONN = atoi(temp);
  my_strncpy(temp,&the_line[25],5);      
  NATTACH = atoi(temp);
  my_strncpy(temp,&the_line[45],5);      
  LABEL = atoi(temp);
  my_strncpy(temp,&the_line[54],1);      
  NDC = atoi(temp);
  for (i = 0; i < LABEL; i++)
    fgets(the_line,sizeof(the_line),file1);
  if (NDC == 3)
    for (i = 1; i <= Atoms; i++)
      fgets(the_line,sizeof(the_line),file1);
  conn = 0;
/*  printf("NCONN = %d NATTCH = %d \n",NCONN,NATTACH); */
  for (i = 0; i < NCONN; i++)
  {
    fgets(the_line,sizeof(the_line),file1);
    tokens = count_tokens(the_line,"\t\n ");
    for (j = 0; j < (tokens - 1); j++)
    {
      a = atoi(gettoken(the_line,"\t\n ",j+1));
      b = atoi(gettoken(the_line,"\t\n ",j+2));
      if ((a != 0) && (b != 0))
      {
	Start(conn) = a;
	End(conn) = b;
	conn++;
      }
    }
  } 
  attach = conn;
  while (attach < conn + NATTACH)
  {
    fgets(the_line,sizeof(the_line),file1);
    tokens = count_tokens(the_line,"\t\n ");
    for (j = 0; j < tokens; j+=2)
    {
      Start(attach) = atoi(gettoken(the_line,"\t\n ",j+1));
      End(attach) = atoi(gettoken(the_line,"\t\n ",j+2));
      attach++;
    }
  } 
  Bonds = attach;
  i = 1;
  while (i <= Atoms)
  {
    UpdateProgress();
    fgets(the_line,sizeof(the_line),file1);
    tokens = count_tokens(the_line,"\t\n ");
    if (tokens > 4)
    {
      sscanf(the_line,"%lf%lf%lf%s",
	     &X(i),&Y(i),&Z(i),Type(i));
      sscanf(&the_line[40],"%lf%lf%lf%s",
	     &X(i+1),&Y(i+1),&Z(i+1),Type(i+1));
      i+=2;
    }
    else
    {
      sscanf(the_line,"%lf%lf%lf%s",
	     &X(i),&Y(i),&Z(i),Type(i));
      i++;
    }
  }  

  column = locate_input_type("MM2");
  for (i = 1; i <= Atoms; i++)
  {
    strcpy(temp,Type(i));
    Atomic_number(i) = get_input_type(i,column,temp,Type(i),dummy);    
  }
  dissect_connection_table(mol);
  assign_bond_order(mol);
  return(TRUE);
}







