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

FILE : wrpdb.c
AUTHOR(S) : Pat Walters
DATE : 1-93
PURPOSE : Routines to write a Brookhaven PDB type file.

******/

static int res_count = 1;

#include "bbltyp.h"

int
  write_dock_pdb(FILE *file1, ums_type *mol)
{
  int i,j;
  char type_name[5], padded_name[5];
  int result;
  char the_res[5], dock_name[5];
  int dock_type;
  int res_num;

  fprintf(file1,"%-51s%-51s\n",Title,Title);

  for (i = 1; i <= Atoms; i ++)
  {
    result =  get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    result =  get_output_type(i,"DOK",Type(i),dock_name,all_caps);
    dock_type = atoi(dock_name);

    strcpy(the_res,"UNK");
    sprintf(padded_name,"%2s",type_name);
    strcpy(type_name,padded_name);
    res_num = res_count;

    fprintf(file1,"ATOM   %4d %-5s%3s%3s%6d%9.3f%8.3f%8.3f%8.3f%8.3f%3d\n",
	    i,
	    type_name,
	    the_res,
	    "",
	    res_num,
	    X(i),
	    Y(i),
	    Z(i),
	    Charge(i),
	    0.0,
	    dock_type);
  }
  fprintf(file1,"TER\n");
  res_count++;
  return(TRUE);
}







