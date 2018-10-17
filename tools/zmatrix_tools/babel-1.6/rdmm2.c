/*****
  This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
  
  For more information please contact :
  
  babel@mercury.aichem.arizona.edu
  ---------------------------------------------------------------------------
  
  FILE : read_mm2.c
  AUTHOR(S) : Pat Walters
  DATE : 10-92
  PURPOSE : Routines to read an mm2 output file
  
  ******/

#include "bbltyp.h"

int 
  read_mm2(FILE *file1, ums_type *mol)
{
  int i,j = 0;
  char the_line[BUFF_SIZE];
  double strain;
  char temp1[5],temp2[5];
  int start, end;
  int is_odd = FALSE;
  int column;
  
  fgets(the_line,sizeof(the_line),file1);
  fgets(the_line,sizeof(the_line),file1);
  sscanf(the_line,"%lf",&strain);
  fgets(the_line,sizeof(the_line),file1);
  sscanf(the_line,"%d%d",&Atoms,&Bonds);

  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  Energy = strain;
  for (i= 1; i <= Atoms; i ++)
    Valence(i) = 0;
  
  if (Atoms % 2 != 0)
  {
    is_odd = TRUE;
    Atoms -= 1;
  }

  column = locate_input_type("MM2");
  for (i= 1; i <= Atoms; i+=2)
  {
    UpdateProgress();
    j = i + 1;
    fgets(the_line,sizeof(the_line),file1);
    sscanf(the_line,"%lf%lf%lf%s%lf%lf%lf%s",
	   &X(i),&Y(i),&Z(i),temp1,&X(j),&Y(j),&Z(j),temp2);

    Atomic_number(i) = get_input_type(i,column,temp1,Type(i),dummy);    
    Atomic_number(j) = get_input_type(j,column,temp2,Type(j),dummy);    
  }

  if (is_odd)
  {
    Atoms += 1;
    fgets(the_line,sizeof(the_line),file1);
    sscanf(the_line,"%lf%lf%lf%s",&X(Atoms),&Y(Atoms),&Z(Atoms),temp2);
    Atomic_number(Atoms) = get_input_type(Atoms,column,temp2,Type(Atoms),dummy);    
  }

  ShowProgress(Bonds,"Reading Bonds");
  j = 0;
  for (i = 0; i < Bonds; i ++)
  {
    UpdateProgress();
    fgets(the_line,sizeof(the_line),file1);
    sscanf(the_line,"%d%d",&start,&end);
    Start(j) = start;
    End(j) = end;
    j++;
  }
  assign_bond_order(mol);
  dissect_connection_table(mol);
  Bonds = j;
  return(TRUE);  
}













