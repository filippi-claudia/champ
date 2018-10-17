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

FILE : wrmm3.c
AUTHOR(S) : Pat Walters
DATE : 2-18-94
PURPOSE : Routines to write mm3 files
******/


#include "bbltyp.h"

int 
write_mm3(FILE *file1, ums_type *mol)
{ 
  int i;
  int type_name;
  int connections = 0;
  int attachments = 0;
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
  double DEL; /* Termianation of geometry optimization */
  int NSYMM;/* Number of symmetry matrices */
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
  METHOD = 1;
  N = Atoms;
  IPRINT = 3;
  if (isdigit(OutputKeywords[0]))
    IPRINT = atoi(OutputKeywords);
  NSTR = 0;
  INIT = 0;
  NCONST = 0;
  TMAX = 999.0;
/*------ CARD 2 -------*/
  DEL = 0.00008;
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
  

  fprintf(file1,"%-60s%d%4d %d  %d %d  %d%-5.0f\n",
	  ID,
	  METHOD,
	  N,
	  IPRINT,
	  NSTR,
	  INIT,
	  NCONST,
	  TMAX);

  fprintf(file1,"%1d%4d%5s%4.5f%8s%5d%5d%5d%5d%5d%5d%5d%5d%10d%5d\n",
	  0,
	  NCON,
	  "",
	  DEL,
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
    get_output_type(i,"MM2",Type(i),temp_type,all_caps);
    type_name = atoi(temp_type);
    type_name = update_mm2_types(mol,i,type_name);
    fprintf(file1,"  %8.5f  %8.5f  %8.5f%5d(%3d)\n",
	    X(i),
	    Y(i),
	    Z(i),
	    type_name,
	    i);
  }
  return(TRUE);
}
    










