/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
------------------------------------------------------------------------------

FILE : wrxyz.c
AUTHOR(S) : Pat Walters
DATE : 11-92
PURPOSE : Routines to write an MSI .BGF file
******/

#include "bbltyp.h"

int 
write_bgf(FILE *file1, ums_type *mol)
{ 
  int i,j;
  char elmnt_typ[5], dreid_typ[5], atm_sym[10], max_val_str[5];
  int max_val;

  assign_hybrid_radii(mol);
  dearomatize(mol);
  
  fprintf(file1,"BIOGRF 200\n");
  fprintf(file1,"DESCRP %s\n",Title);
  fprintf(file1,"REMARK BGF file created by Babel %s\n",BABEL_VERSION);
  fprintf(file1,"FORCEFIELD DREIDING  \n");
  fprintf(file1,"FORMAT ATOM   (a6,1x,i5,1x,a5,1x,a3,1x,a1,1x,a5,3f10.5,1x,a5,i3,i2,1x,f8.5)\n");
  for (i = 1; i <= Atoms; i++)
  {
    get_output_type(i,"XYZ",Type(i),elmnt_typ,all_caps);
    get_output_type(i,"DRE",Type(i),dreid_typ,all_caps);
    get_output_type(i,"HAD",Type(i),max_val_str,all_caps);
    max_val = atoi(max_val_str);
    if (max_val == 0)
      max_val = 1;
    sprintf(atm_sym,"%s%d",elmnt_typ,i);
    fprintf(file1,"%6s %5d %-5s %3s %1s %5s%10.5f%10.5f%10.5f %-5s%3d%2d %8.5f\n",
	    "HETATM",
	    i,
	    atm_sym,
	    "RES",
	    "A",
	    "444",
	    X(i),
	    Y(i),
	    Z(i),
	    dreid_typ,
	    max_val,
	    0,
	    Charge(i));
  }
  fprintf(file1,"FORMAT CONECT (a6,12i6)\n");
  for (i = 1; i <= Atoms; i ++)
  {
    if (Valence(i) > 0)
      {
	fprintf(file1,"CONECT%6d",i);
	for (j = 0; j < Valence(i); j ++)
	  {
	    fprintf(file1,"%6d",Connection(i,j));
	  }
	fprintf(file1,"\n");
	fprintf(file1,"ORDER %6d",i);
	for (j = 0; j < Valence(i); j ++)
	  {
	    fprintf(file1,"%6d",BO(i,j));
	  }
	fprintf(file1,"\n");
      }
  }
  fprintf(file1,"END\n");
  return(TRUE);
}


  
      









