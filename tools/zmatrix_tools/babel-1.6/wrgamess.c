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

FILE : wrgau.c
AUTHOR(S) : Pat Walters
DATE : 1-93
PURPOSE : Routines to write a GAMESS input file

******/



#include "bbltyp.h"

int write_gamess_input(FILE *file1, ums_type *mol)
{
  uppercase(OutputKeywords);
  if (EQ(OutputKeywords,"ZMT"))
      write_gamess_zmatrix(file1,mol);
      else
      if (EQ(OutputKeywords,"CART"))
	  write_gamess_cart(file1,mol);
	  else
	  if (EQ(OutputKeywords,"ZMTMPC"))
	      write_gamess_mopac(file1,mol);
	      else
	      write_gamess_cart(file1,mol);
  return(TRUE);
}

int 
write_gamess_cart(FILE *file1, ums_type *mol)
{ 
  int i;
  char type_name[5];
  int result;
  int at_num;
  
  fprintf(file1," $CONTRL COORD=CART UNITS=ANGS $END\n");
  fprintf(file1,"$DATA\n");
  fprintf(file1,"%s\n",Title);
  fprintf(file1,"Put symmetry info here\n\n");
  
  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    at_num = get_atomic_number(type_name);
    fprintf(file1,"%-3s %4.1f    %8.5f  %8.5f  %8.5f \n",
	    type_name,
	    (float)at_num,
	    X(i),
	    Y(i),
	    Z(i));
  }
  fprintf(file1,"$END\n");
  return(TRUE);
}


int 
write_gamess_zmatrix(FILE *file1, ums_type *mol)
{
  int i=0;
  char type_name[5];
  double a,d;
  int result;
  
  if (mol->internal == NULL)
  {
    initialize_internal(&mol);
    cartint(mol);  
    cartgeom(mol);
  }

  fprintf(file1,"$CONTROL COORD=ZMT NZMAT=0 UNITS=ANGS $END\n");
  fprintf(file1,"$DATA\n");
  fprintf(file1,"%s\n",Title);
  fprintf(file1,"Put symmetry info here\n\n");
  
  for(i=1; i <= Atoms; i++)
  {
    strcpy(type_name,Type(i));
    result = xlate_std_type("XYZ",Type(i),type_name);
    if (result == 0)
    {
      fprintf(stderr,"Unable to assign XYZ type to atom %d type = %s\n",
	      i,Type(i));
      strcpy(type_name,Type(i));
      uppercase(type_name);
    }
    switch(i)
    {
    case 1 :
      fprintf(file1,"  %s \n",type_name);
      break;
    case 2 :
      fprintf(file1,"  %s    %d r%d \n",
	      type_name,
	      mol->internal[i].na,i);
      break;
    case 3 :
      fprintf(file1,"  %s    %d r%d    %d a%d\n",
	      type_name,
	      mol->internal[i].na,i,
	      mol->internal[i].nb,i);
      break;
    default :
      fprintf(file1,"  %s    %d r%d    %d a%d    %d d%d\n",
	      type_name,
	      mol->internal[i].na,i,
	      mol->internal[i].nb,i,
	      mol->internal[i].nc,i); 
    }
  }
  
fprintf(file1,"\n");

  for(i=2; i <= Atoms; i++)
  {
    a = fix_g90_angle(mol->internal[i].w);
    d = fix_g90_angle(mol->internal[i].t);
    switch(i)
    {
    case 2 :
      fprintf(file1,"r2= %6.4f\n",mol->internal[2].r);
      break;
    case 3 :
      fprintf(file1,"r3= %6.4f\na3= %6.2f\n",
	      mol->internal[2].r, a);
      break;
    default :
	fprintf(file1,"r%d= %6.4f\na%d= %6.2f\nd%d= %6.2f\n",
		i,mol->internal[i].r,
		i,a,
		i,d);
      break;
    }
  }
  fprintf(file1,"$END\n");
  return(1);
}

int 
  write_gamess_mopac(FILE *file1, ums_type *mol)
{
  int i=0;
  int num1=0,num2=0,num3=0;
  char type_name[5];
  int default_movements = TRUE;
  int result;
  
  
  if (mol->internal == NULL)
  {
    initialize_internal(&mol);
    cartint(mol);
    cartgeom(mol);
  }

  fprintf(file1,"$CONTROL COORD=ZMTMPC UNITS=ANGS $END\n");
  fprintf(file1,"$DATA\n");
  fprintf(file1,"%s\n",Title);
  fprintf(file1,"Put symmetry info here\n\n");
  
  if (mol->atoms[0].pos[0] == -1)
    default_movements = TRUE;
    
  /* since default atom movements are to be assigned, atoms 1 2 and 3 are
     treated separately */

  for(i=1; i <= Atoms; i++)
  {
    if (default_movements == TRUE)
      switch(i)
      {
      case 1:
	num1=0;
	num2=0;
	num3=0;
	break;
      case 2:
	num1=1;
	num2=0;
	num3=0;
	break;
      case 3:
	num1=1;
	num2=1;
	num3=0;
	break;
      default:
	num1=1;
	num2=1;
	num3=1;
      }
    else
    {
      num1 = mol->atoms[i].pos[0];
      num2 = mol->atoms[i].pos[1];
      num3 = mol->atoms[i].pos[2];
    }
    
    result = xlate_std_type("XYZ",Type(i),type_name);
    if (result == 0)
    {
      fprintf(stderr,"Unable to assign XYZ type to atom %d type = %s\n",
	      i,Type(i));
      strcpy(type_name,Type(i));
      uppercase(type_name);
    }
    fprintf(file1," %-2s%12.7f  %d  %12.6f  %d  %12.6f  %d  %3d  %3d  %3d      %4.4f\n",
	    type_name,
	    mol->internal[i].r,num1,
	    mol->internal[i].w,num2,
	    mol->internal[i].t,num3,
	    mol->internal[i].na,
	    mol->internal[i].nb,
	    mol->internal[i].nc,0.0); 
  }

  fprintf(file1,"$END\n");
  return(TRUE);
}




