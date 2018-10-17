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

FILE : rdxyz.c
AUTHOR(S) : Pat Walters
DATE : 4-94
PURPOSE : routines to read the Quanta format used by MSI
NOTES :  This file format is FORTRAN BINARY.  Each write is surrounded by
the number of bytes being written.  The PRE and POST macros are there to 
"gobble up" this extra information.
	
******/

#include "bbltyp.h"

#ifndef AICHEM 
static warning wstr;
#endif

int dmmy;
#undef DEBUG

#define PRE fread(&dmmy,sizeof(int),1,file1)
#define POST fread(&dmmy,sizeof(int),1,file1)

typedef struct
{
  short int order;
  short int start;
  short int end;
} bnd_type;



int read_quanta(FILE *file1, ums_type *mol)
{

  int i=0,j;
  float the_x,the_y,the_z;
  char the_name[10];
  short int res_num;
  short int atm_typ;
  float charge;
  int nseg, ngroup, natom;
  char version[11];
  char header[201];
  char seg_name[5];
  char title[81];
  short int res_ptr;
  char res_id[7];
  char res_name[5];
  short int atm_ptr;
  short int grp_atm_cnt;
  short int seg_num;
  float xcenter, ycenter, zcenter, rad;
  short int *grp_array;
  char3 *quanta_types;
  int index;
  
  quanta_types = (char3 *)malloc(250 * 3 * sizeof(char));
  if (!read_quanta_types(quanta_types))
  {
    free(quanta_types);
    return(FALSE);
  }

  if (!file1)
  {
    show_warning("Error opening file");
    return(FALSE);
  }
  PRE;
  fread(&nseg,sizeof(int),1,file1);
  fread(&ngroup,sizeof(int),1,file1);
  fread(&natom,sizeof(int),1,file1);
  fread(version,10 * sizeof(char),1,file1);
  fread(header,200 * sizeof(char),1,file1);
  POST;

  grp_array = (short int *)malloc(ngroup * sizeof(short int));

#ifdef DEBUG
  printf("segments = %d\n",nseg);
  printf("groups = %d\n",ngroup);
  printf("atoms = %d\n",natom);
  printf("version = %s\n",version);
  printf("header = %s\n",header);
#endif

  Atoms = natom;
  initialize_ums(&mol);

  strcpy(title,"");
  while (strstr(title,"END") == NULL) 
  {
    PRE; 
    i++;
    fread(title,80 * sizeof(char),1,file1);
    POST;
  }
  PRE;
  for (i = 0; i < nseg; i++)
  {  
    fread(seg_name,4 * sizeof(char),1,file1);
  }
  POST;
  PRE;
  for (i = 0; i < nseg; i++)
  {
    fread(&res_ptr,sizeof(short int),1,file1);
  }
  POST;
  PRE;
  for (i = 0; i < nseg; i++)
  {
    fread(&grp_array[i],sizeof(short int),1,file1);
  }
  POST;
  
  /*  Card 4 - Group Data  */
  PRE;
  for (i = 0; i < nseg; i++)
  {
    for (j = 0; j < grp_array[i]; j++)
    {
      fread(res_id,6 * sizeof(char),1,file1);
    }
  }
  POST;
  PRE;
  for (i = 0; i < nseg; i++)
  {
    for (j = 0; j < grp_array[i]; j++)
    {
      fread(res_name,4 * sizeof(char),1,file1);
    }
  }
  POST;
  PRE;
  for (i = 0; i < nseg; i++)
  {
    for (j = 0; j < grp_array[i]; j++)
    {
      fread(&atm_ptr,sizeof(short int),1,file1);
    }
  }
  POST;
  PRE;
  for (i = 0; i < nseg; i++)
  {
    for (j = 0; j < grp_array[i]; j++)
    {
      fread(&grp_atm_cnt,sizeof(short int),1,file1);
    }
  }
  POST;
  PRE;
  for (i = 0; i < nseg; i++)
  {
    for (j = 0; j < grp_array[i]; j++)
    {
      fread(&seg_num,sizeof(short int),1,file1);
    }
  }
  POST;
  PRE;
  for (i = 0; i < nseg; i++)
  {
    for (j = 0; j < grp_array[i]; j++)
    {
      fread(&xcenter,sizeof(float),1,file1);
    }
  }
  POST;
  PRE;
  for (i = 0; i < nseg; i++)
  {
    for (j = 0; j < grp_array[i]; j++)
    {
      fread(&ycenter,sizeof(float),1,file1);
    }
  }
  POST;
  PRE;
  for (i = 0; i < nseg; i++)
  {
    for (j = 0; j < grp_array[i]; j++)
    {
      fread(&zcenter,sizeof(float),1,file1);
    }
  }
  POST;
  PRE;
  for (i = 0; i < nseg; i++)
  {
    for (j = 0; j < grp_array[i]; j++)
    {
      fread(&rad,sizeof(float),1,file1); 
    }
  }
  POST;

  /* Card 5 - Atom Data */
  PRE;
  for (i = 1; i <= Atoms; i++)
  {
    fread(&the_x,sizeof(float),1,file1);
    X(i) = (double)the_x;
  }
  POST;
  PRE;
  for (i = 1; i <= natom; i++)
  {
    fread(&the_y,sizeof(float),1,file1);
    Y(i) = (double)the_y;
  }
  POST;
  PRE;
  for (i = 1; i <= natom; i++)
  {
    fread(&the_z,sizeof(float),1,file1);
    Z(i) = (double)the_z;
  }
  POST;
  PRE;
  for (i = 1; i <= natom; i++)
  {
    fread(the_name,6 * sizeof(char),1,file1);
    the_name[6] = '\0';
    strcpy(Type(i),the_name);
    clean_atom_type(Type(i));
  }
  POST;
  PRE;
  for (i = 1; i <= natom; i++)
  {
    fread(&res_num,sizeof(short int),1,file1);
  }
  POST;
  PRE;
  for (i = 1; i <= natom; i++)
  {
/*    strcpy(the_name,Type(i));*/
    fread(&atm_typ,sizeof(short int),1,file1);
    index = (int)atm_typ;
/*    strcpy(Type(i),quanta_types[index]);*/
    if EQ(Type(i),"XX")
      printf("Cant assign type %d for atom %d %s\n",index,i,the_name);
  }
  POST;
  PRE;
  for (i = 1; i <= natom; i++)
  {
    fread(&charge,sizeof(float),1,file1);
  }
  POST;
  assign_radii(mol);
  assign_bonds(mol);
  assign_types(mol);
  build_connection_table(mol);
  assign_bond_order(mol);
  free(quanta_types);
  free(grp_array);
  read_to_eof(file1);
  return(TRUE);
}

int read_quanta_types(char3 *quanta_types)
{
  FILE *file1;
  char the_line[80];
  int i;
  char the_type[3];
  char babel_dir[80];

  for (i = 0; i < 250; i++)
    strcpy(quanta_types[i],"XX");
  
  if ((file1 = fopen("quanta.lis","r")) == NULL) 
  {
    if (getenv("BABEL_DIR") == NULL)
    {
      show_warning("The environment variable BABEL_DIR is not defined");
      show_warning("Please define this variable to so Babel can find quanta.lis\n");
      return(FALSE);
    }
    else 
      strcpy(babel_dir,getenv("BABEL_DIR"));
    strcat(babel_dir,"/quanta.lis");  
    if ((file1 = fopen(babel_dir,"r")) == NULL) 
    {
      sprintf(wstr,"Could not open Quanta types files file %s \n",babel_dir);
      show_warning(wstr);
      return(FALSE);
    }
  }
  while (fgets(the_line,sizeof(the_line),file1))
  {
    sscanf(the_line,"%d %s",&i,the_type);
    strcpy(quanta_types[i],the_type);
  }
}

  
    
	 
	   
  


