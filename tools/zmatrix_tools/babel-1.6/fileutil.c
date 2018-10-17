/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
-------------------------------------------------------------------------------

FILE : file_util.c
AUTHOR(S) : Pat Walters
DATE : 11-95
PURPOSE : routines for dealing with files
******/

#include "bbltyp.h"

static char wstr[BUFF_SIZE];

FILE *open_w_env(char *f_name, char *env_var)
{
  FILE *fp;
  char new_name[BUFF_SIZE];
  char message[BUFF_SIZE];
  
  fp = fopen(f_name,"r");
  if (fp)
    return(fp);
  else
     if (getenv(env_var))
     {
        strcpy(new_name,getenv(env_var));
        strcat(new_name,"/");
        strcat(new_name,f_name);
        if ((fp = fopen(new_name,"r")) != NULL)
           return(fp);
     }

  sprintf(message,"Unable to open %s - Please define the %s environment variable",
          f_name,env_var);
  fatal_error(message);
  
  return(fp);
}


FILE *open_read(char *file_name)
{
  FILE *file1 = NULL;

  if ((file1 = fopen(file_name,"r")) == NULL)
  {
    sprintf(wstr,"unable to open file %s\n",file_name);
    fatal_error(wstr);
    exit(0);
  }
  return(file1);
}

FILE *open_write(char *file_name)
{
  FILE *file1 = NULL;
  
  if EQ(file_name,"CON")
    file1 = stdout;
  else
    if ((file1 = fopen(file_name,"w")) == NULL)
    {
      printf("unable to open file %s as output\n",file_name);
      exit(0);
    }
  return(file1);
}

void close_file(char *file_name, FILE *file1)
{
  if NOTEQn(file_name,"CON",3)
    fclose(file1);
}



