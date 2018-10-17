/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : rdmdl.c
AUTHOR(S) : Pat Walters
DATE : 1-6-94
PURPOSE : routines to read an MDL mol file

******/

#include "bbltyp.h"
#define DELIMS "\t\n "

int 
read_mdl(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i,result,found;
  char temp[15];
  char temp_title[100];
  long int pos,prev;
  double dummyd;
  char dummyc[50];
  int dummyi;

  strcpy(temp_title,InfileName);

#ifdef OLD_MDL_READER
  fgets(the_line,sizeof(the_line),file1);
  if (count_tokens(the_line,DELIMS) > 0)
    strcpy(temp_title,gettoken(the_line,DELIMS,1));
  else
    strcpy(temp_title,InfileName);
  for (i = 0; i < 2; i++)
  {
    fgets(the_line,sizeof(the_line),file1);
  }
  fgets(the_line,sizeof(the_line),file1);
  my_strncpy(temp,the_line,3);
  Atoms = atoi(temp);
  my_strncpy(temp,&the_line[3],3);
  Bonds = atoi(temp);
  ShowProgress(Atoms,"Reading Atoms");
  result = initialize_ums(&mol);
  strcpy(Title,temp_title);
#endif

  prev = ftell(file1);
   do
   {
      fgets(the_line,sizeof(the_line),file1);
      if  (sscanf(the_line,"%lf%lf%lf%s%ld%ld%ld%ld%ld",
                  &dummyd,&dummyd,&dummyd,dummyc,
                  &dummyi,&dummyi,&dummyi,&dummyi,&dummyi,&dummyi) == 9)
         if (isalpha(dummyc[0]))
            break;

      prev = pos;
      pos = ftell(file1);
   } while (!check_for_eof(file1));

  if (check_for_eof(file1))
  {
     show_warning("Unable to find beginning of file");
     return(FALSE);
  }
  
  fseek(file1,prev,0);
  fgets(the_line,sizeof(the_line),file1);
  sscanf(&the_line[3],"%3d",&Bonds);
  the_line[3] = '\0';
  sscanf(the_line,"%3d",&Atoms);

  ShowProgress(Atoms,"Reading Atoms");
  result = initialize_ums(&mol);
  if (!result)
  {
    show_warning("Unable to initialize ums");
    Atoms = Bonds = 0;
    return(FALSE);
  }

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
  dissect_connection_table(mol);
  assign_types(mol);

#ifdef OLD_MDL_READER

  pos = ftell(file1);
  fgets(the_line,sizeof(the_line),file1);
  found = FALSE;
  if ((strstr(the_line,"END")) || (strstr(the_line,"$$$")))
    found = TRUE;
  if (!found)
    fseek(file1,pos,0);

#endif

   prev = ftell(file1);
   while (fgets(the_line,sizeof(the_line),file1))
   {
      if  (sscanf(the_line,"%lf%lf%lf%s%d%d%d%d%d%d",
                  &dummyd,&dummyd,&dummyd,dummyc,
                  &dummyi,&dummyi,&dummyi,&dummyi,&dummyi,&dummyi) == 9)
         if (isalpha(dummyc[0]))
            break;
   } 

  if (!check_for_eof(file1))
     fseek(file1,prev,0);


  return(TRUE);
}



double get_scale_factor(ums_type *mol)
{
  int i;
  double dist;
  
  for (i = 0; i < Bonds; i++)
    if ((Type(Start(i))[0] == 'C') && (Type(End(i))[0] == 'C'))
      if (Bond_order(i) == 1)
    {
      dist = distance(Point(Start(i)),Point(End(i)));
      return(dist);
    }
  return(0.0);
}


void scale_for_ChemWindow(ums_type *mol)
{
  double fact;
  int i;
  
  fact = get_scale_factor(mol);
  if (fact == 0.0)
  {
    show_warning("Couldn't find a C - C single bond to use for scaling");
    show_warning("The bond lengths are probably going to be screwed up");
  }
  else
    for (i = 1; i <= Atoms; i++)
    {
      X(i) = X(i) * 1.54/fact;
      Y(i) = Y(i) * 1.54/fact;
      Z(i) = Z(i) * 1.54/fact;
    }
}


  

 



   
    
    
    
	  


