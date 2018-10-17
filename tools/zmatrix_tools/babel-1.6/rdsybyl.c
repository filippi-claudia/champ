/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
------------------------------------------------------------------------

FILE : rdsybyl.c
AUTHOR(S) : Jussi Eloranta <eloranta@tukki.jyu.fi>
	    University of Jyvdskyld, Finland

Modified:  6-94 by Matthew Stahl to work with new method of 
           handling substructure files

Modified:  8-95 by Pat Walters to accomodate the modified mol2
           format used by the NCI Molecules R' Us folks

DATE : 12-93
PURPOSE : routines to read an SYBYL MOL2 file
          Supported record types (@<tripos>):
	  molecule, atom, bond. Others are igonred.

******/

#include "bbltyp.h"

#define MOL "MOLECULE"
#define ATO "ATOM"
#define BON "BOND"

int
read_sybyl(FILE *file1, ums_type *mol)
{
  char input_line[BUFF_SIZE];
  char temp_title[BUFF_SIZE];
  int i;
  char temp_type[5];
  int start,end,order;
  int found_string = FALSE;
  long int pos;
  char bo_str[5];
  double temp_num = 0.0;
  int column;

  while (fgets(input_line,sizeof(input_line),file1)) 
  {
    if ((input_line[0] == '@') && strstr(input_line,MOL))
    {
      found_string = TRUE;
      break;
    }
  }

  if (!found_string)
    return(FALSE);

  /* @<tripos>molecule record (6 data lines) */

  /* molecule name */
  fgets(input_line,sizeof(input_line),file1);
/*
  if(sscanf(input_line, "%s",temp_title) != 1)	
    show_warning("Problem reading name of current file");
*/
  strcpy(temp_title,input_line);

  /* # atoms, # bonds, # substructures, # features, # sets */
  fgets(input_line,sizeof(input_line),file1);
  if(sscanf(input_line, "%d %d", &Atoms, &Bonds) != 2)
  {
    show_warning("Problem reading number of Atoms and Bonds");
    return(FALSE);
  }

  found_string = FALSE;
  while (fgets(input_line, sizeof(input_line), file1))
  {
    if (strstr(input_line,"Energy"))
    {
      sscanf(input_line,"%*s%*s%lf",&temp_num);
    }

    if ((input_line[0] == '@') && strstr(input_line,ATO))
    {
      found_string = TRUE;
      break;
    }
  }
  
  if (!found_string)
  {
    show_warning("Unable to locate Atom Table in Sybyl MOL2 file");
    Atoms = Bonds = 0;
    return(FALSE);
  }
  
  ShowProgress(Atoms,"Reading Atoms");
  
  initialize_ums(&mol);
  initialize_residues(&mol);
  strcpy(Title,temp_title);
  strip_return(Title);
  Energy = temp_num;

  column = locate_input_type("SYB");  
  for (i = 1; i <= Atoms; i ++) 
  {
    UpdateProgress();

    strcpy(ResName(i),"UNK");
    ResNum(i) = 1;
    ChainNum(i) = 0;
    fgets(input_line,sizeof(input_line),file1);
    sscanf(input_line," %*s %s %lf %lf %lf %s %*s %*s %lf",
	   AtmId(i),
	   &X(i),
	   &Y(i),
	   &Z(i),
	   temp_type,
	   &Charge(i));
    Atomic_number(i) = get_input_type(i,column,temp_type,Type(i),dummy);    
  }
  
  found_string = FALSE;
  while (fgets(input_line, sizeof(input_line), file1))
  {
    if ((input_line[0] == '@') && strstr(input_line,BON))
    {
      found_string = TRUE;
      break;
    }
  }

  if (!found_string)
  {
    show_warning("Unable to locate Bond Table in Sybyl MOL2 file");
    release_ums(mol);
    free(mol->residues);
    Atoms = Bonds = 0;
    return(FALSE);
  }
  
  for (i = 0; i < Bonds; i ++) 
    {
      fgets(input_line,sizeof(input_line),file1);
      sscanf(input_line,"%*d %d %d %s",&start,&end,bo_str);
      uppercase(bo_str);
      if (EQ(bo_str,"AR"))
	order = 5;
      else 
	if (EQ(bo_str,"AM"))
	  order = 1;
	else
	  order = atoi(bo_str);
      Connection(start,Valence(start)) = end;
      BO(start,Valence(start)) = order;
      Valence(start)++;      
      Connection(end,Valence(end)) = start;
      BO(end,Valence(end)) = order;
      Valence(end)++;
    }

  build_connection_table(mol);

  pos = ftell(file1);

  while (fgets(input_line,sizeof(input_line), file1) != NULL) 
  {
    if ((input_line[0] == '@') && strstr(input_line,MOL))
    {
      fseek(file1,pos,0);
      break;
    }
  }

  type_hydrogens(mol);

  return(TRUE);
}




