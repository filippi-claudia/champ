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

FILE : rdisis.c
AUTHOR(S) : Pat Walters
DATE : 1-6-94
PURPOSE : routines to read an MDL isis file

******/

#include "bbltyp.h"

int 
read_isis(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE] = {'\0'};
  char temp_title[BUFF_SIZE] = {'\0'};
  char previous[BUFF_SIZE] = {'\0'};
  int i = 0;
  int result = FALSE;
  char temp[15] = {'\0'};
  int temp_num;

  while (!check_for_eof(file1)) 
  {
    fgets(the_line,sizeof(the_line),file1);
    my_strncpy(temp,&the_line[20],2);
    if (EQn(temp,"3D",2) || EQn(temp,"2D",2))
    {
      strcpy(temp_title,the_line);
      strip_return(temp_title);
      break;
    }
    strcpy(previous,the_line);
  }
  fgets(the_line,sizeof(the_line),file1);
/* End of changes ... Ajay, Nov. 17 1994! */

  fgets(the_line,sizeof(the_line),file1);
  my_strncpy(temp,the_line,3);
  Atoms = atoi(temp);
  my_strncpy(temp,&the_line[3],3);
  Bonds = atoi(temp);
  
  ShowProgress(Atoms,"Reading Atoms");
  result = initialize_ums(&mol);
  strip_return(previous);
  strcpy(Title,previous);
  for (i = MIN_ATOM; i <= Atoms; i++)
    {
      UpdateProgress();
      fgets(the_line,sizeof(the_line),file1);
      my_strncpy(temp,the_line,10);
      X(i) = atof(temp);
      my_strncpy(temp,&the_line[10],10);
      Y(i) = atof(temp);
      my_strncpy(temp,&the_line[20],10);
      Z(i) = atof(temp);
      my_strncpy(Type(i),&the_line[31],2);
      clean_atom_type(Type(i));
    }
  ShowProgress(Atoms,"Reading Atoms");
  for (i = 0; i < Bonds; i++)
    {
      UpdateProgress();
      fgets(the_line,sizeof(the_line),file1);
      my_strncpy(temp,the_line,3);
      Start(i) = atoi(temp);
      my_strncpy(temp,&the_line[3],3);
      End(i) = atoi(temp);
      my_strncpy(temp,&the_line[6],3);
      Bond_order(i) = atoi(temp);
      if (Bond_order(i) == 4)
	Bond_order(i) = 5;
    }
/* 
  The following line has been added to finish reading the molecule.
  That is, read till "$$$$" line. Ajay ... Nov. 17, 1994
*/
    while (!check_for_eof(file1)) 
    {
      fgets(the_line,sizeof(the_line),file1);
      if (EQn(the_line,">  <MOLREGNO>",13))
      {
	fgets(the_line,sizeof(the_line),file1);
	sscanf(the_line,"%d",&temp_num);
	sprintf(Title,"MOL_%03d",temp_num);
      }
      
      if ((EQn(the_line,">  <MDLNUMBER>",14)) || 
	  (EQn(the_line,">  <mdlnumber>",14)) || 
	  (EQn(the_line,">  <VXNUMBER>",13)) || 
	  (EQn(the_line,">  <MODELREGNO>",14)))
      {
	fgets(the_line,sizeof(the_line),file1);
	strip_return(the_line);
	my_strncpy(Title,the_line,25);
      }
      my_strncpy(temp,the_line,4);
      if (strncmp(temp,"$$$$",4) == 0)
	break;
    }

  dissect_connection_table(mol);
  
  assign_type_by_bo(mol);
  return(TRUE);
}  

 



   
    
    
    
	  


