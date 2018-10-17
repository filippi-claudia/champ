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

FILE : rdxed.c
AUTHOR(S) : Keith Trollope (kit1000@cam.ac.uk)
DATE : 2-95
PURPOSE :  Read a XED (new COSMIC) file
******/
#include "bbltyp.h"

int 
read_xed(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i, ptr;
  int line_count=1;
  char temp_type[10],mass_str[10];
  int multi_file = FALSE;
  int mass;
  float chge;
  int column;
  ums_type *new_mol;  /* necessary in case have to delete XEDS or `fields' */

  fgets(the_line,sizeof(the_line),file1);
  sscanf(the_line,"%10lf%10d%10d",&Energy,&Atoms,&Bonds);
  ShowProgress(Atoms,"Reading Atoms");
  initialize_ums(&mol);
  fgets(the_line,sizeof(the_line),file1);
  sscanf(the_line,"%s",Title);
  for(i=0;i<Bonds;i++)
  {
    if (!(i%5)) 
    {
      fgets(the_line,sizeof(the_line),file1);
      ptr=0;
    }
    sscanf(&the_line[ptr],"%d%d",&Start(i),&End(i));
    ptr += 16;
  }

  dissect_connection_table(mol);
  column = locate_input_type("XED");
  for(i=1;i<=Atoms;i++)
  {
    fgets(the_line,sizeof(the_line),file1);
    sscanf(the_line,"%6s%15lf%15lf%15lf%6s%12lf",
	   mass_str,
	   &X(i),
	   &Y(i),
	   &Z(i),
	   temp_type,
	   &chge);
    mass=atoi(mass_str);
    if (mass < 1) 
      strcpy(temp_type,"21"); /* want to remove these funny points */
    Atomic_number(i) = get_input_type(i,column,temp_type,Type(i),dummy);
  }  
  fgets(the_line,sizeof(the_line),file1);
  new_mol = delete_atoms(mol,"X");
  *mol = *new_mol;
  return(TRUE);
}


void xed_add_connection(ums_type *mol, int start, int end, int bnd_num)
{
    Connection(start,Valence(start))=end;
    Valence(start)++;
    Connection(end,Valence(end))=start;
    Valence(end)++;
    Start(bnd_num)=start;
    End(bnd_num)=end;
}




