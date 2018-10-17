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

FILE : wrsybyl.c
AUTHOR(S) : Pat Walters based on original version by
            Jussi Eloranta <eloranta@tukki.jyu.fi>
	    University of Jyvdskyld, Finland
DATE : 12-93
PURPOSE : Write a SYBYL mol2 file
******/

#include "bbltyp.h"

int write_sybyl(FILE *file1, ums_type *mol)
{
  int i;
  char temp_type[5];
  char ele[5];
  char bond_ord[3];
  pdb_type_rec *pdb_types = NULL;
  char type_label[10];
  char res_label[10];
  char res_num[10];

  if (HasResidues == FALSE)
    pdb_types = (pdb_type_rec *)malloc((Atoms + 1) * sizeof(pdb_type_rec));


  uppercase(OutputKeywords);
  if (NOTEQ(OutputKeywords,"NOAROM"))
    find_aromatic_atoms(mol);

  fprintf(file1, "@<TRIPOS>MOLECULE\n");
  if (strlen(Title) == 0)
    strcpy(Title,"*****");
  fprintf(file1, "%s\n",Title);
  fprintf(file1, " %d %d 0 0 0\n", Atoms, Bonds);
  fprintf(file1, "SMALL\n");
  fprintf(file1, "NO_CHARGES\n\n");
/*
  fprintf(file1, "Energy = %10.3f\n",Energy);
*/
  fprintf(file1, "@<TRIPOS>ATOM\n");
  for(i = 1; i <= Atoms; i++) {
    get_output_type(i,"SYB",Type(i),temp_type,dummy);
    get_element_type(mol,i,ele);

    if (HasResidues)
    {
      strcpy(type_label,trim_spaces(AtmId(i)));
/*      strcpy(res_label,trim_spaces(ResName(i))); */
      strcpy(res_label,"<1>");
      sprintf(res_num,"%d",ResNum(i));
    }
    else
    {
      strcpy(pdb_types[i].name,ele);
      assign_pdb_number(pdb_types,i);
      sprintf(type_label,"%s%d",ele,pdb_types[i].number);
      strcpy(res_label,"<1>");
      strcpy(res_num,"1");
    }
    
    fprintf(file1,"%7d%1s%-6s%12.4f%10.4f%10.4f%1s%-5s%4s%1s %-8s%10.4f\n",
	    i,"",type_label, X(i), Y(i), Z(i),"",temp_type,res_num,"",res_label,Charge(i));
  }

  fprintf(file1, "@<TRIPOS>BOND\n");
  for(i = 0; i < Bonds; i++)
  {
    if (Bond_order(i) != 5)
      sprintf(bond_ord,"%d",Bond_order(i));
    else
      strcpy(bond_ord,"ar");
    fprintf(file1, "%6d%6d%6d%3s%2s\n", i+1, Start(i), End(i),"",bond_ord);
  }
  if (pdb_types)
    free(pdb_types);
  fprintf(file1,"\n");
  return(TRUE);
}


void fix_sybyl_types(ums_type *mol)
{
  int i,j,conn;
  int is_OCo2, is_C_cat;
  
  for (i = 1; i <= Atoms; i++)
  {
    switch(Atomic_number(i))
    {
    case 8 :
      if (EQ(Type(i),"O-"))
      {
	is_OCo2 = FALSE;
	for (j = 0; j < Valence(i); j++)
	{
	  conn = Connection(i,j);
	  /* In Sybyl atoms should be typed as O.co2 only if they
	     are in carboxylate or phosphate groups */
	  
	  if (EQ(Type(conn),"Cac") || EQ(Type(conn),"Sac") || EQ(Type(conn),"Pac"))
	  {
	    is_OCo2 = TRUE;
	    break;
	  }
	}
	if (!is_OCo2)
	  strcpy(Type(i),"O2");
      }
      break;
    }
  }
}









