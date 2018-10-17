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

FILE : read_types.c
AUTHOR(S) : Pat Walters
DATE : 1-93
PURPOSE : routines to read the types table

******/

#include "bbltyp.h"

static warning wstr;

static char ***type_table = NULL;
static int rows;
static int columns;

#define SIZE 8

void read_types_table()
{
  FILE *file1;
  int i,j;
#ifndef AICHEM
  char babel_dir[80];
#endif
  char the_line[BUFF_SIZE];
  int tokens;

  if ((file1 = fopen("types.lis","r")) == NULL) 
  {
#ifndef AICHEM 
    if (getenv("BABEL_DIR") == NULL)
    {
      show_warning("The environment variable BABEL_DIR is not defined");
      show_warning("Please define this variable to so babel can find element.lis\n");
      exit(0);
    }
    else 
      strcpy(babel_dir,getenv("BABEL_DIR"));
    strcat(babel_dir,"/types.lis");  
    if ((file1 = fopen(babel_dir,"r")) == NULL) 
    {
      sprintf(wstr,"Could not open types files file %s \n",babel_dir);
      fatal_error(wstr);
    }
#else
    file1 = open_read("/usr/local/babel/types.lis");
#endif 
  }
  fscanf(file1,"%d %d",&rows,&columns);
  
  type_table = (char ***) malloc (rows * sizeof (char *) );
  if (type_table) 
  {
    for (i = 0; i < rows; i++) 
    {
      type_table[i] = (char **) malloc (columns * sizeof(char*) );		
    }
  }
  
  for (i = 0; i < rows; i++)
    for (j = 0; j < columns; j++)
    {
      type_table[i][j] = (char *) malloc(SIZE * sizeof(char));
    }
  
  fgets(the_line,sizeof(the_line),file1);  
  for (i = 0; i < rows; i++)
  {
  fgets(the_line,sizeof(the_line),file1);  
    tokens = count_tokens(the_line,"\t\n ");
    if (tokens != columns)
    {
      show_warning("There is a problem with the file types.lis");
      sprintf(wstr,"The problem appears to be with line %d",i+2);
      show_warning(wstr);
      show_warning(the_line);
      sprintf(wstr,"This line should have %d tokens, but it has %d tokens",
	      columns,tokens);
      show_warning(wstr);
      exit(0);
    }

    for (j = 0; j < columns; j++)
      get_token(type_table[i][j],the_line,"\n\t ",(j+1));
  }
  fclose(file1);
}

void clean_up_types_table(void)
{
  int i,j;
  
  for (i = 0; i < rows; i++)
    for (j = 0; j < columns; j++)
    {
      free(type_table[i][j]);
    }
  for (i = 0; i < rows; i++) 
    {
      free(type_table[i]);
    }
  free(type_table);
}


void write_type_table(void)
{
  int i,j;
  
  for (i = 0; i < rows; i++)
  {
    for (j = 0; j < columns; j++)
      printf("%5s",type_table[i][j]);
    printf("\n");
  }
}

/*-----------------------------------------
FUNCTION : locate_input_type

look for a format string in row 0 of type_table
return column number where format string was
found
-------------------------------------------*/

int locate_input_type(char *format)
{
  int col_num;
  int i;

  if (type_table == NULL)
  {
    sprintf(wstr,"No type table -- you forgot to call babel_init()");
    fatal_error(wstr);
  }
  
  for (i = 0; i < columns; i++)
    if EQ(type_table[0][i],format)
      col_num = i;
  if (col_num == -1)
  {
    sprintf(wstr,"Unknown format %s\n",format);
    fatal_error(wstr);
  }
  return(col_num);
}

int get_input_type(int atom_num,int col_num, char *input, char *std_type, enum type_err error)
{
  int i;
  int atomic_num;
  int found = FALSE;
  
  for (i = 1; i < rows; i++)
    if EQ(type_table[i][col_num],input)
    {
      strcpy(std_type,type_table[i][0]);
      atomic_num = atoi(type_table[i][1]);
      found = TRUE;
      break;
    }
  
  if (!found)
  {
    sprintf(wstr,"Unable to assign UMS type to atom %d type = %s\n",
	    atom_num,input);
    show_warning(wstr);
    switch(error)
    {
    case zero :
      strcpy(std_type,"0");
      break;
    case dummy :
      strcpy(std_type,"Du");
      break;
    case all_caps :
      uppercase(input);
      strcpy(std_type,input);
      break;
    default :
      strcpy(std_type,input);
      break;
    }
  }
  return(atomic_num);
}


int get_output_type(int at_num,char *format, char *input, char *out_type, enum type_err error)
{
  int i;
  int col_num = -1;
  int found = FALSE;

  for (i = 0; i < columns; i++)
    if EQ(type_table[0][i],format)
      col_num = i;
  if (col_num == -1)
  {
    sprintf(wstr,"Unknown format %s\n",format);
    fatal_error(wstr);
  }

  for (i = 1; i < rows; i++)
    if EQ(type_table[i][0],input)
    {
      strcpy(out_type,type_table[i][col_num]);
      found = TRUE;
      break;
    }

  if (!found)
  {
    sprintf(wstr,"Unable to assign %s type to atom %d type = %s",
	    format,at_num,input);
    show_warning(wstr);
    switch(error)
    {
    case zero :
      strcpy(out_type,"0");
      break;
    case dummy :
      get_output_type(i,format,"X",out_type,all_caps);
      break;
    case all_caps :
      uppercase(input);
      strcpy(out_type,input);
      break;
    default :
      strcpy(out_type,input);
      break;
    }
  }
  return(found);
}

int get_std_type(char *format, char *input, char *std_type)
{
  int i;
  int col_num = -1;
  
  for (i = 0; i < columns; i++)
    if EQ(type_table[0][i],format)
      col_num = i;
  if (col_num == -1)
  {
    sprintf(wstr,"Unknown format %s\n",format);
    fatal_error(wstr);
  }
  for (i = 1; i < rows; i++)
    if EQ(type_table[i][col_num],input)
    {
      strcpy(input,type_table[i][0]);
      return(1);
    }
  return(0);
}

int xlate_std_type(char *format, char *std_type, char *output)
{
  int i;
  int col_num = -1;
  
  for (i = 0; i < columns; i++)
    if EQ(type_table[0][i],format)
      col_num = i;
  if (col_num == -1)
  {
    sprintf(wstr,"Unknown format %s\n",format);
    fatal_error(wstr);
  }
  for (i = 1; i < rows; i++)
    if EQ(type_table[i][0],std_type)
    {
      strcpy(output,type_table[i][col_num]);
      return(1);
    }
  return(0);
}

  


