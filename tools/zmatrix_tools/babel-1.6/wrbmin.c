/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
------------------------------------------------------------------------------

FILE : wrbmin.c
AUTHOR(S) : Pat Walters
DATE : 11-94
PURPOSE : write a batchmin com file to minimize a multi-structure file
******/

#include "bbltyp.h"

int 
write_bmin_com(FILE *file1, ums_type *mol)
{ 
  int i;
  char temp_name[BUFF_SIZE];
  int done = FALSE;
  int *heavies;
  int heavy_count = 0;
  
  heavies = (int *)malloc(Atoms * sizeof(int));
  for (i = 1; i <= Atoms; i++)
  {
    if (Type(i)[0] != 'H')
    {
      heavies[heavy_count] = i;
      heavy_count++;
    }
  }    
  fprintf(file1,"%s\n",InfileName);
  strcpy(temp_name,InfileName);
  strip_extension(temp_name,temp_name);
  strcat(temp_name,".reduce");
  fprintf(file1,"%s\n",temp_name);
  fprintf(file1," DEMX       0      0      0      0    50.0000     0.0000");
  for(i = 0;i < heavy_count; i++)
  {
    if (((i % 4) == 0) && (i < heavy_count ))
    {
      fprintf(file1,"\n COMP ");
      done = TRUE;
    }
    fprintf(file1,"%7d",heavies[i]);    
  }
  fprintf(file1,"\n");
  fprintf(file1," MULT\n");
  fprintf(file1," FFLD       1\n");
  fprintf(file1," BGIN \n");
  fprintf(file1," READ\n");
  fprintf(file1," MINI       1      0  10000\n");
  fprintf(file1," END\n");
  return(TRUE);
}


  
      









