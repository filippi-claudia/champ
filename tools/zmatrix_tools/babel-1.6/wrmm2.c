/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
----------------------------------------------------------------------------

FILE : wrmm2.c
AUTHOR(S) : Pat Walters
DATE : 11-92
PURPOSE : Routines to write mm2 input and mm2 output type files, also 
contains a really wierd file format used by Pat's MOUSE program.

******/


#include "bbltyp.h"

int 
write_mm2(FILE *file1, ums_type *mol)
{ 
  int i;
  int type_name;
  char temp_type[5];
  
  fprintf(file1,"%s\n",Title);
  fprintf(file1,"  %8.5f\n",Energy);
  fprintf(file1,"%5d%5d \n",Atoms,Bonds);

  for(i = 1;i <= Atoms; i++)
  {
    get_output_type(i,"MM2",Type(i),temp_type,dummy);
    type_name = atoi(temp_type);
    type_name = update_mm2_types(mol,i,type_name);

    fprintf(file1,"  %8.4f  %8.4f  %8.4f%5d     ",
	    X(i),
	    Y(i),
	    Z(i),
	    type_name);
    if (((i % 2) == 0)||(i == Atoms)) 
      fprintf(file1,"\n");
  }
  
  for(i = 0;i < Bonds; i++)
  {
    fprintf(file1,"%5d%5d\n",
	    Start(i),
	    End(i));
  }
  return(TRUE);
}

int 
write_mouse(FILE *file1, ums_type *mol)
{ 
  int i;
  
  fprintf(file1,"%s\n",Title);
  fprintf(file1,"  %8.5f\n",0.00);
  fprintf(file1,"%5d%5d \n",Atoms,Bonds);
  
  for(i = 1;i <= Atoms; i++)
  {
    fprintf(file1,"  %8.4f  %8.4f  %8.4f%5s     ",
	    X(i),
	    Y(i),
	    Z(i),
	    Type(i));
    if (((i % 2) == 0)||(i == Atoms)) 
      fprintf(file1,"\n");
  }
  
  for(i = 0;i < Bonds; i++)
  {
    fprintf(file1,"%5d%5d\n",
	    Start(i),
	    End(i));
  }
  return(TRUE);
}


int 
write_mm2_input(FILE *file1, ums_type *mol)
{ 
  int i;
  int type_name;
  int connections = 0;
  int attachments = 0;
  int result;
  char temp_type[5];


  char ID[60]; /* filename */
  int METHOD; /* 0 no cojugated pi system, 1 if conjugated pi system */
  int N; /* #of atoms */
  int IPRINT; /* Controls amount of printout */
  int NSTR; /* Restricted motion data  */
  int INIT; /* Minimize energy  */
  int NCONST; /* Read in new constants ? */
  double TMAX; /* Max time */

  int NCON; /* Number of connected atoms */
  int NATTACH; /*Number of attached atoms */
  int NSYMM;/* Number of symmetry matrics */
  int NX; /* Number of coordiante calcualtions or replacement cards */
  int NROT; /* Reorient */
  int LABEL; /* Change names or atomic weights */
  int NDC; /* Dipole and charge interaction energy */
  int NCALC; /* Crystal conversions */
  int HFORM; /* Heat of formation */
  int MVDW; /* Approximate van der Waals */
  int NDRIVE; /* Dihedral driver */
  
  
  for (i = 0;i < Bonds; i++)
  {
    if ((Valence(Start(i)) == 1) || (Valence(End(i)) == 1))
      attachments ++;
  }
  
  connections = Bonds - attachments;  
  strcpy(ID,OutfileName);

/*------ CARD 1 -------*/
  METHOD = 0;
  N = Atoms;
  IPRINT = 3;
  if (isdigit(OutputKeywords[0]))
    IPRINT = atoi(OutputKeywords);
  NSTR = 0;
  INIT = 0;
  NCONST = 0;
  TMAX = 999.0;
/*------ CARD 2 -------*/
  NCON = connections;
  NATTACH = attachments;
  NSYMM = 0;
  NX = 0;
  NROT = 0;
  LABEL = 0;
  NDC = 0;
  NCALC = 0;
  HFORM = 0;
  MVDW = 1;
  NDRIVE = 0;
  

#ifdef AICHEM
  fprintf(file1,"    1\n"); /* This card only necessary for MM2 78 */
#endif
  fprintf(file1,"%-60s%d%4d %d  %d %d  %d%-5.0f\n",
	  ID,
	  METHOD,
	  N,
	  IPRINT,
	  NSTR,
	  INIT,
	  NCONST,
	  TMAX);
  fprintf(file1,"%5d%20s%5d%5d%5d%5d%5d%5d%5d%5d%10d%5d\n",
	   NCON,
	   "",
	   NATTACH,
	   NSYMM,
	   NX,
	   NROT,
	   LABEL,
	   NDC,
	   NCALC,
	   HFORM,
	   MVDW,
	   NDRIVE);

  for(i = 0;i < Bonds; i++)
  {
    if ((Valence(Start(i)) > 1) && (Valence(End(i)) > 1))
      fprintf(file1,"%5d%5d\n",
	    Start(i),
	    End(i));
  }
  
  attachments = 0;
  
  for(i = 0;i < Bonds; i++)
  {
    if ((Valence(Start(i)) == 1) || (Valence(End(i)) == 1))
    {
      attachments ++;
      fprintf(file1,"%5d%5d",
	    Start(i),
	    End(i));

      if (((attachments % 8) == 0))
	fprintf(file1,"\n");
    }
  }
  
  if (((attachments % 8) != 0))
    fprintf(file1,"\n");
  
  for (i = 1;i <= Atoms; i++)
  {
    result = xlate_std_type("MM2",Type(i),temp_type);
    if (result == 0)
    {
      fprintf(stderr,"Unable to assign MM2 type to atom %d type = %s\n",
	      i,temp_type);
      strcpy(Type(i),"0");
    }
    type_name = atoi(temp_type);
    type_name = update_mm2_types(mol,i,type_name);
    fprintf(file1,"  %8.5f  %8.5f  %8.5f%5d     ",
	    X(i),
	    Y(i),
	    Z(i),
	    type_name);
    if (((i % 2) == 0)||(i == Atoms)) 
      fprintf(file1,"\n");
  }
  return(TRUE);
}


int update_mm2_types(ums_type *mol, int i, int temp_type)
{
  int type_name;

  type_name = temp_type;
  switch (temp_type)
  {
  case 5 :
    type_name = type_mm2_hydrogens(mol,i);
    break;
  case 2 :
    type_name = check_for_carbonyl(mol,i);
    break;
  case 1:
    type_name = check_for_cyclopropane(mol,i);
    break;
  }
  return(type_name);
}


int 
check_for_carbonyl(ums_type *mol, int atom_num)
{
  int bonded_atom;
  int i;
  
  for (i = 0; i < Valence(atom_num); i ++)
  {
    bonded_atom = Connection(atom_num,i);
    if (EQ(Type(bonded_atom),"O2") || EQ(Type(bonded_atom),"S2"))
      return(3);
  }
  return(2);
}

int check_for_cyclopropane(ums_type *mol, int atom1)
{
  int i,j;
  
  for (i = 0; i < Valence(atom1); i++)
    for (j = 0; j < Valence(atom1); j++)
    {
      if (bonded(mol,Connection(atom1,i),Connection(atom1,j)))
	return(22);
    }
  return(1);
}


int bonded(ums_type *mol, int atom1, int atom2)
{
  int i;
  
  for (i = 0; i < Valence(atom1); i++)
    if (Connection(atom1,i) == atom2)
      return(TRUE);
  return(FALSE);
}

    

int type_mm2_hydrogens(ums_type *mol, int i)
{
  int j,k,l,result = 0;
  
  j = Connection(i,0);
  if (Type(j)[0] == 'O')
  {
    for (k = 0; k < Valence(j); k++)
    {
      l = Connection(j,k);
      if ((EQ(Type(l),"C2")) || (EQ(Type(l),"Cac")))
	result = 24;
    }
    if (result == 0)
      result = 21;
  }
  if (Type(j)[0] == 'N')
    result = 23;
  if  (result == 0)
    result = 5;
  return(result);
}

    
