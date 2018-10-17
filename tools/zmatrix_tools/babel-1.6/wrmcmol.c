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

FILE : wrmcmol.c
AUTHOR(S) : Matt Stahl
DATE : 7-93
PURPOSE : Routines to write a MacMolecule file 

******/

#include "bbltyp.h"

static warning wstr;

extern element_type *elements; 

int 
write_mcmol(FILE *file1, ums_type *mol)
{ 
  int i,j,num_types = 0;
  int result,unique;
  double vdw_radius,bs_radius,red,grn,blu;  
  char type_name[5];
  char type_list[80];
  char token[5];
  pdb_type_rec *pdb_types;

  strcpy(OutputKeywords,"BS");
  fprintf(file1,";%s\n",Title);
  fprintf(file1,";Stick color\n");
  fprintf(file1,"2(0.700, 0.700, 0.700)\n");
  fprintf(file1,";Model\n");
  
  if (strcmp(OutputKeywords,"KEYWORDS GO HERE") == 0)
  {
    strcpy(OutputKeywords,"BS");
    sprintf(wstr,"Model type not specified: (default is ball and stick)\n");
    show_warning(wstr);
  }
  fprintf(file1,"%s\n",OutputKeywords);

  strcpy(type_list,"");
  for (i = 1; i <= Atoms; i ++)
  {
    unique = TRUE;
    result = get_output_type(i,"MCML",Type(i),type_name,dummy);    
    for (j = 1;j <= num_types;j++)
    {
      strcpy(token,gettoken(type_list,"/",j));
      if (EQ(type_name,token))
      {
	unique = FALSE;
      }
    }

    if (unique)
    {
      num_types++;
      strcat(type_list,type_name);
      strcat(type_list,"/");
      get_atom_info(Type(i),&vdw_radius,&bs_radius,&red,&grn,&blu);
      
      fprintf(file1,"%s=%5.3f,%5.3f(%5.3f,%5.3f,%5.3f)\n",
	      type_name,
	      vdw_radius,
	      bs_radius,
	      red,
	      grn,
	      blu);
    }
  }
  
  pdb_types = (pdb_type_rec *)malloc((Atoms + 1)* 
				     sizeof(pdb_type_rec));

  for (i = 1; i <= Atoms; i ++)
  {
    result = xlate_std_type("MCML",Type(i),type_name);
    if (result == 0)
    {
      strcpy(type_name,"X");
    }

    strcpy(pdb_types[i].name,type_name);
    assign_pdb_number(pdb_types,i);
    fprintf(file1,"%1s%d:  %7.4f  %7.4f  %7.4f\n",
	    type_name,
	    pdb_types[i].number,
	    X(i),
	    Y(i),
	    Z(i));
  }

  fprintf(file1,";Bonds\n");
  for(i = 0; i < Bonds; i++)
  {
    fprintf(file1,"%s%d,%s%d\n",
	    pdb_types[Start(i)].name,
	    pdb_types[Start(i)].number,
	    pdb_types[End(i)].name,
	    pdb_types[End(i)].number);
    
  }
  return(TRUE);
}



void get_atom_info(char *type_name,double *vdw_radius,double *bs_radius, double *red, double *grn, double *blu)
{
  int i;
  int found;
  char trans_name[5];
 
  xlate_std_type("XYZ",type_name,trans_name);
  found = FALSE;

  for (i = 0; i < MAX_ELEMENTS; i++)
  { 
    if (strncmp(trans_name,elements[i].name,2) == 0) 
    {
      *vdw_radius = elements[i].vdw_rad;
      *bs_radius = elements[i].bs_rad;
      *red = elements[i].red;
	  *grn = elements[i].green;
	  *blu = elements[i].blue;
      found = TRUE;
    }
 
  }

  if (found == FALSE)
    {    
     *vdw_radius = 1.5;
     *bs_radius = 0.75;
     *red = 0.00;
     *grn = 0.00;
     *blu = 0.00;
   }
  
}




void translate_color(int *color,double *red,double *grn,double *blu)
{
  switch(*color)
  {
  case BBL_UNDEF:
    *red = 0.8;*grn = 0.0;*blu = 0.6;
    break;
  case BBL_BLACK:
    *red = 0.0;*grn = 0.0;*blu = 0.0;
    break;
  case BBL_GREY:
    *red = 0.5;*grn = 0.5;*blu = 0.5;
    break;
  case BBL_DKBLU:
    *red = 0.05;*grn = 0.05;*blu = 0.5;
    break;
  case BBL_BLUE:
    *red = 0.0;*grn = 0.0;*blu = 1.0;
    break;
  case BBL_LTBLU:
    *red = 0.1;*grn = 0.4;*blu = 1.0;
    break;
  case BBL_AQUA:
    *red = 0.07;*grn = 0.5;*blu = 0.7;
    break;
  case BBL_TURQ:
    *red = 0.07;*grn = 0.7;*blu = 0.7;
    break;
  case BBL_BLUGRN:
    *red = 0.32;*grn = 0.87;*blu = 0.67;
    break;
  case BBL_DKGRN:
    *red = 0.04;*grn = 0.94;*blu = 0.04;
    break;
  case BBL_GREEN:
    *red = 0.0;*grn = 0.0;*blu = 0.1;
    break;
  case BBL_LTGRN:
    *red = 0.4;*grn = 1.0;*blu = 0.1;
    break;
  case BBL_YELGRN:
    *red = 0.8;*grn = 1.0;*blu = 0.1;
    break;
  case BBL_YELLOW:
    *red = 1.0;*grn = 1.0;*blu = 0.0;
    break;
  case BBL_ORANGE:
    *red = 1.0;*grn = 0.37;*blu = 0.08;
    break;
  case BBL_DKRED:
    *red = 0.7;*grn = 0.0;*blu = 0.1;
    break;
  case BBL_RED:
    *red = 1.0;*grn = 0.0;*blu = 0.0;
    break;
  case BBL_PINK:
    *red = 1.0;*grn = 0.1;*blu = 0.5;
    break;
  case BBL_REDPUR:
    *red = 1.0;*grn = 0.1;*blu = 0.8;
    break;
  case BBL_PURPLE:
    *red = 0.75;*grn = 0.08;*blu = 0.5;
    break;
  case BBL_BLUPUR:
    *red = 0.5;*grn = 0.1;*blu = 0.8;
    break;
  case BBL_WHITE:
    *red = 1.0;*grn = 1.0;*blu = 1.0;
    break;
  }
}
