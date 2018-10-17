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

FILE : menus.c
AUTHOR(S) : Pat Walters
DATE : 4-93 (complete rewrite)
PURPOSE : menu interface to babel

******/

#include "bbltyp.h"
#include "bblmast.h"

void babel_menus(ums_type *mol)
{
  Verbose = TRUE;
  read_menu(mol);
  write_menu(mol);
  continuation_menu(mol);
}

void read_menu(ums_type *mol)
{
  int i,j = 0,k;
  int choice;
  int has_keywords = FALSE;

  printf("%s Babel %s %s\n",STARS,BABEL_VERSION,STARS);
  puts("Input file type \n");
  for (i = 0; i < Size; i++)
  {
    if (master[i].operation == input)
    {
      j++;
      printf("%4d. %-20s \t",j,master[i].type_name);
      if (j % 2 == 0)
	printf("\n");
    }
  }
  printf("\n");
  choice = get_choice(1,j);
  k = 0;
  for (i = 0; i < Size; i++)
  {
    if (master[i].operation == input)
    {
      k++;
      if (k == choice)
      {
	InputInfo = master[i];
	has_keywords = TRUE;
      }
    }
  }
  printf("Input file name : ");
  gets(InfileName);
  printf("Keywords : ");
  gets(InputKeywords);

  if (strlen(InputKeywords) < 1)
    strcpy(InputKeywords,"KEYWORDS GO HERE");
  do_inputs(mol);
}


void write_menu(ums_type *mol)
{
  int i,j = 0,k;
  int choice;
  int has_keywords = FALSE;
  FILE *outfile;

  printf("%s Babel %s %s\n",STARS,BABEL_VERSION,STARS);
  puts("Output file type \n");
  for (i = 0; i < Size; i++)
  {
    if (master[i].operation == output)
    {
      j++;
      printf("%4d. %-20s \t",j,master[i].type_name);
      if (j % 2 == 0)
	printf("\n");
    }
  }
  printf("\n");
  choice = get_choice(1,j);
  k = 0;
  for (i = 0; i < Size; i++)
  {
    if (master[i].operation == output)
    {
      k++;
      if (k == choice)
      {
	OutputInfo = master[i];
	has_keywords = TRUE;
      }
    }
  }
  printf("Output file name : ");
  gets(OutfileName);
  if (has_keywords)
  {
    printf("Keywords : ");
    gets(OutputKeywords);
    if (strlen(OutputKeywords) < 1)
      strcpy(OutputKeywords,"KEYWORDS GO HERE");
  }
  outfile = open_write(OutfileName);
  do_outputs(outfile,mol);
  close_file(OutfileName,outfile);
}


int 
continuation_menu(ums_type *mol)
{
  int choice;
  int result;
  int done = FALSE;

  while (done == FALSE)
  {
    puts("\n***************************************");
    puts("1. read a file ");
    puts("2. write a file ");
    puts("3. quit");  
  
    choice = get_choice(1,3);

    switch(choice)
    {
    case 1 :
      result = release_ums(mol);
      read_menu(mol);
      write_menu(mol);
      break;
    case 2 :
      write_menu(mol);
      break;
    case 3 :
      done = TRUE;
      break;
    }
  }
  return(TRUE);
}


int get_choice(int min, int max)
{
  char choice_string[100];
  int choice = -1;
  int done = FALSE;
  int i;
  
  while (done == FALSE)
  {
    printf("Choice : ");
    gets(choice_string);
    for (i = 0; i < (int) strlen(choice_string); i++)
    {
      if (!isdigit(choice_string[i]))
	choice = 0;
    }
    if (choice == -1)
      choice = atoi(&choice_string[0]);
    if ((choice >= min) && (choice <= max))
      return(choice);
    printf("%d is not a valid choice\n",choice);
    choice = -1;
    printf("Please try again\n");
  }
  return(TRUE);
}

      

