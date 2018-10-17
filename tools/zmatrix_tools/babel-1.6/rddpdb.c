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

FILE : rdpdb.c
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : Routines to read a Brookhave Protien Data Bank file

******/
#include "bbltyp.h"



int 
read_dock_pdb(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  char temp_title[BUFF_SIZE];
  char temp_type[10];
  int count = 0;
  int i;
  int result;
  long int pos;
  int column;
 
  pos = ftell(file1);
  while ((fgets(the_line,sizeof(the_line), file1) != NULL) && NOTEQn(the_line,"TER",3))
  {
    if (EQn(the_line,"ATOM",4))
      count++;
  }
  Atoms = count;
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  fseek(file1,pos,0);
  i = MIN_ATOM;

  fgets(the_line,sizeof(the_line), file1);
  sscanf(&the_line[6],"%s",temp_title);

  column = locate_input_type("DOK");  
  while ((fgets(the_line,sizeof(the_line), file1) != NULL) && NOTEQn(the_line,"TER",3))
  {
    if ((EQn(the_line,"REMARK",6)) &&(strstr(the_line,"score")))
    {
      sscanf(the_line,"%*s %*s %*s %*s %lf",&Energy);
    }
	  
    sprintf(Title,"%s %10.3f",temp_title,Energy);

    if (EQn(the_line,"ATOM",4))
    {
      UpdateProgress();
      sscanf(the_line,"%*s%*s %*s %*s%*s %lf%lf%lf%lf%*s%s",
	     &X(i),
	     &Y(i),
	     &Z(i),
	     &Charge(i),
	     temp_type);
      Atomic_number(i) = get_input_type(i,column,temp_type,Type(i),dummy);    
      i++;
    }
  }

  add_element_types(mol);

  if (Atoms > 0)
  {
    result = assign_radii(mol);
    result = assign_bonds(mol);
    result = assign_types(mol);
    assign_bond_order(mol);
  }

  return(TRUE);  
}
