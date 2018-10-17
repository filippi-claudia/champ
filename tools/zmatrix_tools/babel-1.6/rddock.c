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

FILE : rddock.c
AUTHOR(S) : Pat Walters
DATE : 2-29-96 (Leap Day !)
PURPOSE : Routine to read a dock database
	
MODIFIED : 10-16-93 to allow the use of multistructure files

******/


#include "bbltyp.h"

#define DOCK_CHG_FACTOR 1000.0

int first_one = TRUE;

int 
read_dock_database(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  char buffer[BUFF_SIZE];
  int i;
  int column;

  if (first_one)
  {
    fgets(the_line,sizeof(the_line),file1);
    first_one = FALSE;
  }
  
  fgets(the_line,sizeof(the_line),file1);
  my_strncpy(buffer,&the_line[0],50);
  
  fgets(the_line,sizeof(the_line),file1);
  sscanf(the_line,"%d",&Atoms);

  ShowProgress(Atoms,"Reading Atoms");

  initialize_ums(&mol);
  strcpy(Title,buffer);
  rtrim(Title);

  column = locate_input_type("DOK");

  for (i = 1; i <= Atoms; i++)
  {
    UpdateProgress();
    fgets(the_line,sizeof(the_line), file1);

    my_strncpy(buffer,&the_line[0],5);
    X(i) = my_atof(buffer)/1000.0;

    my_strncpy(buffer,&the_line[5],5);
    Y(i) = my_atof(buffer)/1000.0;

    my_strncpy(buffer,&the_line[10],5);
    Z(i) = my_atof(buffer)/1000.0;

    my_strncpy(buffer,&the_line[15],2);

    if (isspace(buffer[0]))
    {
      buffer[0] = buffer[1];
      buffer[1] = '\0';
    }
    Atomic_number(i) = get_input_type(i,column,buffer,Type(i),dummy);

    my_strncpy(buffer,&the_line[17],5);
    Charge(i) = my_atof(buffer)/DOCK_CHG_FACTOR;
  }
  add_element_types(mol);
  assign_radii(mol);
  assign_bonds(mol); 
  assign_types(mol);
  build_connection_table(mol); 
  assign_bond_order(mol);
  return(TRUE);
}

#undef DOCK_CHG_FACTOR
