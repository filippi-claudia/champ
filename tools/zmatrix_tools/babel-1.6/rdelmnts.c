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

FILE : rdelmnts.c
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : routines to read the element list

Modified 4-93 to use the BABEL_DIR environment variable
******/

#include "bbltyp.h"

element_type *elements;

int 
read_element_file()
{
  FILE *file1;
  char element_line[BUFF_SIZE];
  int i;
  char babel_dir[80];
  element_type *temp;

  temp  = (element_type *)malloc(MAX_ELEMENTS * sizeof(element_type));
  elements = temp;
  
  if ((file1 = fopen("element.lis","r")) == NULL) 
  {
#ifndef AICHEM 
    if (getenv("BABEL_DIR") == NULL)
    {
      printf("The environment variable BABEL_DIR is not defined\n");
      printf("Please define this variable to so babel can find element.lis\n");
      exit(0);
    }
    else 
      strcpy(babel_dir,getenv("BABEL_DIR"));
    strcat(babel_dir,"/element.lis");
    if ((file1 = fopen(babel_dir,"r")) == NULL) 
#else
      if ((file1 = fopen("/usr/local/babel/element.lis","r")) == NULL) 
#endif
      {
	printf("Could not open element file %s \n",babel_dir);
	return(0);
      }
  }
  
  while (fgets(element_line,sizeof(element_line), file1) != NULL)
  {
    if (count_tokens(element_line,"\t\n ") != 10)
    {
      printf("Error reading element file\n");
      return(0);
    }
    if (sscanf(element_line,"%d",&i) !=1)
    {
      printf("Error reading element file\n");
      return(0);
    }
    elements[i].number = i;
    sscanf(element_line,"%*d%s%lf%lf%lf%lf%d%lf%lf%lf",
	   elements[i].name,
	   &elements[i].cov_rad,
	   &elements[i].bond_ord_rad,
	   &elements[i].vdw_rad,
	   &elements[i].bs_rad,
	   &elements[i].max_bonds,
	   &elements[i].red,
	   &elements[i].green,
	   &elements[i].blue);
  }
  fclose(file1);
  return(1);
}


void clean_up_elements()
{
  free(elements);
}
    
int 
write_elements(element_type elements[])
{
  int i;
  
  for (i = 0; i < MAX_ELEMENTS; i++)
  {
    printf("%2d [%c][%c] %8.3f %8.3f %8.3f %5d %d\n",
	   elements[i].number,
	   elements[i].name[0],
	   elements[i].name[1],
	   elements[i].bond_ord_rad,
	   elements[i].cov_rad,
	   elements[i].vdw_rad,
	   elements[i].max_bonds,
	   elements[i].color);

  }
  return(TRUE);
}


int get_max_bonds(int atomic_number)
{
  return(elements[atomic_number].max_bonds);
}

void atomic_number_to_name(int i, char *name)
{
  strcpy(name,elements[i].name);
}


/*-----------------------------------------

       FUNCTION - add_element_types

take a ums which only has atomic numbers and
add the element symbols.
----------------------------------------*/

void add_element_types(ums_type *mol)
{
  int i;
  
  for (i = 1; i <= Atoms; i++)
  {
    strcpy(Type(i),elements[Atomic_number(i)].name);
  }
}















