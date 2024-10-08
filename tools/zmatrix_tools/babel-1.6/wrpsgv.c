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

FILE : wrpsgv.c
AUTHOR(S) : Pat Walters, Tom Pollard
DATE : 1-96
PURPOSE : Routines to write a PS-GVB Cartesian file
******/

#include "bbltyp.h"

int
write_psgvb_cart(FILE *file1, ums_type *mol)
{
  int i;
  char type_name[5];
  int result;
  FILE *file2;
  char buffer[BUFF_SIZE];
  char babel_dir[BUFF_SIZE];


  if ((file2 = fopen("psgvb.hdr","r")) != NULL)
  {
    while (fgets(buffer,sizeof(buffer),file2))
      fprintf(file1,"%s",buffer);

    if (file2)
      fclose(file2);
  }
  else
    if (getenv("BABEL_DIR"))
    {
      strcpy(babel_dir,getenv("BABEL_DIR"));
      strcat(babel_dir,"/psgvb.hdr");

      if ((file2 = fopen(babel_dir,"r")) != NULL)
      {
        while (fgets(buffer,sizeof(buffer),file2))
          fprintf(file1,"%s",buffer);

        if (file2)
          fclose(file2);
      }
    }
    else
    {
      fprintf(file1,"inv0230\n");
      fprintf(file1,"JOB: JOBNAME\n");
      fprintf(file1,"DIR: EXEDIR TEMPDIR\n");
    }

    fprintf(file1,"This PS-GVB input file generated by Babel %s\n",BABEL_VERSION);
    fprintf(file1,"$zmat\n");

  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    clean_atom_type(type_name);
    fprintf(file1,"%s%d  %12.7f  %12.7f  %12.7f\n",
       type_name,i,X(i),Y(i),Z(i));
  }

  if (strcmp(OutputKeywords,"KEYWORDS GO HERE") == 0)
    strcpy(OutputKeywords," ");
  fprintf(file1,"$\n$gen\n%s\n$\n",OutputKeywords);

  return(TRUE);
}
