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

FILE : rdmacmol.c
AUTHOR(S) : Matthew Stahl
DATE : 8-93
PURPOSE : routines to read the Mac Molecule format

******/
#include "bbltyp.h"
       
int 
read_mcmol(FILE *file1, ums_type *mol)
{
  char mcmol_line[BUFF_SIZE];
  char temp_type[10];
  char start[5],end[5];
  char from[5],to[5];
  int i,j,count=0;
  int result;
  int l;
  
  while(fgets(mcmol_line,sizeof(mcmol_line),file1) != NULL)
  {
    strcpy(mcmol_line,clean_comments(mcmol_line));
    if (strpbrk(mcmol_line,":") != NULL)
      count++;
  }
  rewind(file1);
  
  Atoms = count;
  ShowProgress(Atoms,"Reading Atoms");
  result = initialize_ums(&mol);

  { 
    i = 1;
    while (fgets(mcmol_line,sizeof(mcmol_line),file1) != NULL)
    {
      UpdateProgress();
      strcpy(mcmol_line,clean_comments(mcmol_line));
      if (strpbrk(mcmol_line,":") != NULL)
      {
	strcpy(Type(i),gettoken(mcmol_line,": \n\t",1));
	X(i) = (double)atof(gettoken(mcmol_line,": \n\t",2));
	Y(i) = (double)atof(gettoken(mcmol_line,": \n\t",3));
	Z(i) = (double)atof(gettoken(mcmol_line,": \n\t",4));
	i++;
      }
    }
    rewind(file1);
    
    i = 0;
    while (fgets(mcmol_line,sizeof(mcmol_line),file1) != NULL)
    {
      strcpy(mcmol_line,clean_comments(mcmol_line));
      if (strpbrk(mcmol_line,",") != NULL &&
	  strpbrk(mcmol_line,"(") == NULL &&
	  strpbrk(mcmol_line,")") == NULL)
      {
	strcpy(start,gettoken(mcmol_line,", \n\t",1));
	strcpy(end,gettoken(mcmol_line,", \n\t",2));
	
	for (j = 1;j <= Atoms;j++)
	{
	 l = strlen(Type(j));
	 if (EQn(Type(j),start,l))
	   Start(i) = j;
	 if (EQn(Type(j),end,l))
	   End(i) = j;
       }
	i++;
      }
      Bonds = i;
    }
  }
  
  for (i = 1;i <= Atoms;i++)
  {
    Type(i)[1] = NULL_CHAR;
  }
  
  if (NOTEQ(InputKeywords,NOKEY))
  {
    printf("Input Keywords = %s\n",InputKeywords);
    if (strpbrk(InputKeywords,"/") == NULL)
      show_warning("Format for specifying atom identities: \"macmolecule atom/output atom\"");
    else
    {
      count = count_tokens(InputKeywords," ");
      for (i = 1;i <= count;i++)
      {
	strcpy(temp_type,gettoken(InputKeywords," ",i));
	strcpy(from,gettoken(temp_type,"/",1));
	strcpy(to,gettoken(temp_type,"/",2));
	for (j = 1; j <= Atoms;j++)
	{
	  if (EQ(Type(j),from))
	    strcpy(Type(j),to);
	}
      }
    }  
  }
      
  dissect_connection_table(mol);
  result = assign_types(mol);
  result = build_connection_table(mol);
  assign_bond_order(mol);
  return(TRUE);
}

char *
  clean_comments(char *mcmol_str)
{
  int i = 0;
  
  while (mcmol_str[i] != NULL_CHAR)
  {
    if (mcmol_str[i] == ';')
    {
      mcmol_str[i] = NULL_CHAR;
      break;
    }
    i++;
  }
  
  return(mcmol_str);
}







   
    
    
    
	  


