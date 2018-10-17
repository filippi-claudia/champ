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

FILE : rdmopac.c
AUTHOR(S) : Pat Walters
DATE : 1-93
PURPOSE : Routines to read a MOPAC output file

******/

#include "bbltyp.h"

int 
read_mopac_output(FILE *file1, ums_type *mol)
{
  char mop_line[BUFF_SIZE];
  int i = 0;
  int j;
  int result;
  long int pos = 0;
  long int chg_pos = 0;
  double theEnergy;
  int has_charges = FALSE;
  int count = 0;
  char trigger1[300], trigger2[300];
  int foundEnergy = FALSE;

  while (fgets(mop_line,sizeof(mop_line),file1))
  {
    if (strstr(mop_line,"FINAL") && strstr(mop_line,"HEAT"))
    {
      sscanf(mop_line,"%*s%*s%*s%*s%*s%lf",&theEnergy);
      foundEnergy = TRUE;
    }

    if (strstr(mop_line,"NET ATOMIC CHARGES"))
    {
      has_charges = TRUE;
      chg_pos = ftell(file1);
    }

    if (strstr(mop_line,"CARTESIAN COORDINATES") && foundEnergy)
      break;
  }

  for (i = 0;i < 3;i++)
    fgets(mop_line,sizeof(mop_line),file1);
  
  pos = ftell(file1);
  
  while (fgets(mop_line,sizeof(mop_line), file1))
  {
    j = sscanf(mop_line,"%s %s",trigger1,trigger2);
    if ((count_tokens(mop_line," ") == 5) && isdigit(trigger1[0]) &&
	isalpha(trigger2[0]))
      count++;
    else
      break;
  }
  fseek(file1,pos,0);

  Atoms = count;
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  Energy = theEnergy;

  for (i = 1; i <= Atoms; i++)
  {
    UpdateProgress();
    fgets(mop_line,sizeof(mop_line), file1);
    j = sscanf(mop_line,"%*d %s %lf %lf %lf",Type(i),&X(i),&Y(i),&Z(i)); 
    if (j != 4)
    {
      show_warning("Input file error");
      release_ums(mol);
      Atoms = 0;
      return(FALSE);
    }
  }
  
  pos = ftell(file1);

  if (has_charges)
  {
    fseek(file1,chg_pos,0);
    toss(file1,2);
    for (i = 1; i <= Atoms; i++)
    {
      fgets(mop_line,sizeof(mop_line), file1);
      sscanf(mop_line,"%*s%*s%lf",&Charge(i));
    }
  }

  if (Atoms > 0)
  {
    result = assign_radii(mol);
    result = assign_bonds(mol);
    result = assign_types(mol);
    result = build_connection_table(mol);
    assign_bond_order(mol);
  }

  fseek(file1,pos,0);

  while (fgets(mop_line,sizeof(mop_line), file1))
  {
    pos = ftell(file1);
    if (strstr(mop_line,"FINAL") && strstr(mop_line,"HEAT"))
    {
      fseek(file1,pos,0);
      break;
    }
  }  
  
  for (i = 1; i <= Atoms; i++)
  {
    fgets(mop_line,sizeof(mop_line), file1);
    sscanf(mop_line,"%*s%*s%f",&Charge(i));
  }
  
  return(TRUE);
}
