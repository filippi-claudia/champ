/****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact:

babel@mercury.aichem.arizona.edu
-------------------------------------------------------------------------------
FILE : rdpcmod.c
AUTHOR(S) : Abby Parrill
DATE : 6-94
PURPOSE : Read a pcmodel file into the UMS

*****/

#include "bbltyp.h"
#include "bblmacs.h"

int read_pcmodel(FILE *file1, ums_type *mol)
{
  char pcmod_line[BUFF_SIZE];
  int i;        /* dummy variable for undesired sections in record */
  char temp[BUFF_SIZE];    /* temporary storage for pcmodel atom type */
  int column;
  
  while (fgets(pcmod_line, sizeof(pcmod_line), file1) != NULL)
  {
    if (strstr(pcmod_line,"NA"))
      break;
  }
  sscanf(pcmod_line, "%*s %d", &Atoms);
  initialize_ums(&mol);    

  column = locate_input_type("PCM");
  while (fgets(pcmod_line, sizeof(pcmod_line), file1) != NULL)
  {
    if (EQn(pcmod_line, "AT", 2))        /* find atom records */
    {
      sscanf(pcmod_line,"%*s %d ",&i);
      get_token(temp,pcmod_line,",: \n\t",3);
      Atomic_number(i) = get_input_type(i,column,temp,Type(i),dummy);    
      get_token(temp,pcmod_line,",: \n\t",4);
      X(i) = atof(temp);
      get_token(temp,pcmod_line,",: \n\t",5);
      Y(i) = atof(temp);
      get_token(temp,pcmod_line,",: \n\t",6);
      Z(i) = atof(temp);
      get_pcmod_bonds(pcmod_line,mol,i);
    }
    if (strchr(pcmod_line,'}'))
    {
      break;
    }
  }
  build_connection_table(mol);
  
  return(TRUE);
}

void get_pcmod_bonds(char *the_line, ums_type *mol, int i)
{
  char *start, *old_start;
  int done = FALSE;
  char delims[] = ", \t\n";
  int tokens,j,k;
  char temp[BUFF_SIZE];
  
  start = strchr(the_line,'B');
  old_start = start;
  old_start++;
  while ((start != '\0') && (!done))
  {
    switch(*start)
    {
    case 'S'  :
    case 'C'  :
    case '\n' :
      *start = '\0';
      done = TRUE;
      break;
    default   :
      start++;
    }
  }
  tokens = count_tokens(old_start,delims);
  Valence(i) = tokens/2;
  k = 0;
  for (j = 1; j < (2 * Valence(i)); j+=2)
  {
    get_token(temp,old_start,delims,j);
    Connection(i,k) = atoi(temp);
    get_token(temp,old_start,delims,j+1);
    BO(i,k) = atoi(temp);
    k++;
  }
}
  




