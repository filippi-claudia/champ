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

FILE : report.c
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : routines to write out the statistics of a molecule (torsions, distances, angles ...)

******/



#include "bbltyp.h"

int print_report_file(FILE *file1, ums_type *mol)
{
  fprintf(file1,"FILENAME: %s\n",InfileName);
  fprintf(file1,"INTERATOMIC DISTANCES \n");
  distance_matrix(mol,file1);
  fprintf(file1,"\n\nBOND ANGLES\n");
  print_angles(mol,file1);
  fprintf(file1,"\n\nTORSION ANGLES\n");
  print_torsions(mol,file1);
  return(TRUE);
}

void get_element_type(ums_type *mol, int i, char *type)
{
  int result;
  char type_name[5];
  
  result = xlate_std_type("XYZ",Type(i),type_name);
    if (result == 0)
    {
      fprintf(stderr,"Unable to assign XYZ type to atom %d type = %s\n",
	      i,Type(i));
      strcpy(type_name,Type(i));
      toupper(type_name[0]);
      tolower(type_name[1]);
    }
  strcpy(type,type_name);
}

void distance_matrix(ums_type *mol,FILE *file1)
{
  int columns = 7;
  int max, min = 1;
  int i,j;
  char type[5];

  max = columns;
  while (max <= Atoms + columns)
  {
    fprintf(file1,"\n");
    get_element_type(mol,min,type);
    fprintf(file1,"%15s%4d",type,min);
    for (i = min + 1; ((i < max) && (i < Atoms)); i++)
      if (i <= Atoms)
      {
	get_element_type(mol,i,type);
	fprintf(file1,"%7s%4d",type,i);
      }
    fprintf(file1,"\n");

    fprintf(file1,"%14s","");
    for (i = min; i < max; i++)
      if (i <= Atoms)
	fprintf(file1,"%11s","-----------");
      
    fprintf(file1,"\n");
    for (i = min; i <= Atoms; i++)
    {
      get_element_type(mol,i,type);
      fprintf(file1,"%4s%4d ",type,i);
      for (j = min; j < max; j++)
	if (j <= i)
	  fprintf(file1,"%10.3f ",distance(Point(i),Point(j)));
      fprintf(file1,"\n");
    }
    max += columns - 1;
    min += columns - 1;
  }
  fprintf(file1,"\n");
}


void print_torsions(ums_type *mol,FILE *file1)
{
  int a,b,c,d;
  int i,j,k;
  double angle1;
  int angle_count = 0;
  torsion_rec *tr;


  tr = (torsion_rec *)malloc(Atoms * 10 * sizeof(torsion_rec));
  if (tr == NULL)
  {
    printf("Memory Allocation Error\n");
    exit(0);
  }
  
  for (i = 0; i < Bonds; i++)
  {
    b = Start(i);
    c = End(i);
    for (j = 0; j < Valence(Start(i)); j ++)
      if (Connection(Start(i),j) != End(i))
      {
	a = Connection(Start(i),j);
	for (k = 0; k < Valence(End(i)); k ++)
	  if ((Connection(End(i),k) != Start(i)) &&
	      (Connection(End(i),k) != a))
	  {
	    d = Connection(End(i),k);
	    tr[angle_count].a = a;
	    tr[angle_count].b = b;
	    tr[angle_count].c = c;
	    tr[angle_count].d = d;
	    angle_count ++;
	  }
      }
  }
  qsort(tr,angle_count,sizeof(torsion_rec),QSORT_PROTO compare_torsion);
  for (i = 0; i < angle_count; i++)
  {
    angle1 = torsion(Point(tr[i].a),Point(tr[i].b),Point(tr[i].c),Point(tr[i].d));
    fprintf(file1,"%4d %4d %4d %4d %10.3f\n",tr[i].a,tr[i].b,tr[i].c,tr[i].d,angle1);
  }
  
  free(tr); 
}	
  


void print_angles(ums_type *mol,FILE *file1)
{
  int a,b,c,d;
  int i,j,k;
  double angle1;
  int angle_count = 0;
  angle_rec *tr;
  angle_rec last;

  tr = (angle_rec *)malloc(Atoms * 10 * sizeof(angle_rec));
  if (tr == NULL)
  {
    printf("Memory Allocation Error\n");
    exit(0);
  }

  last.a = -1;
  last.b = -1;
  last.c = -1;
  
  for (i = 0; i < Bonds; i++)
  {
    b = Start(i);
    c = End(i);
    for (j = 0; j < Valence(Start(i)); j ++)
      if (Connection(Start(i),j) != End(i))
      {
	a = Connection(Start(i),j);
	for (k = 0; k < Valence(End(i)); k ++)
	  if (Connection(End(i),k) != Start(i))
	  {
	    d = Connection(End(i),k);
	    tr[angle_count].a = a;
	    tr[angle_count].b = b;
	    tr[angle_count].c = c;
	    sort_values(&tr[angle_count].a,&tr[angle_count].c);
	    angle_count ++;
	    tr[angle_count].a = b;
	    tr[angle_count].b = c;
	    tr[angle_count].c = d;
	    sort_values(&tr[angle_count].a,&tr[angle_count].c);
	    angle_count ++;
	  }
      }
  }
  qsort(tr,angle_count,sizeof(angle_rec),QSORT_PROTO compare_angle); 
  for (i = 0; i < angle_count; i++)
  {
    if ((tr[i].a != last.a) || (tr[i].b != last.b) || (tr[i].c != last.c)) 
    {
      angle1 = bond_angle(Point(tr[i].a),Point(tr[i].b),Point(tr[i].c));
      fprintf(file1,"%4d %4d %4d %4s %4s %4s %10.3f \n",
	      tr[i].a,tr[i].b,tr[i].c,
	      Type(tr[i].a),Type(tr[i].b),Type(tr[i].c),angle1);
      last.a = tr[i].a;
      last.b = tr[i].b;
      last.c = tr[i].c;
    }
  }
  free(tr);
}


void sort_values(int *a, int *b)
{
  int temp;
  
  if (*a > *b)
  {
    temp = *a;
    *a = *b;
    *b = temp;
  }
}

  

int compare_torsion(torsion_rec *t1, torsion_rec *t2)
{
  if (t1->a > t2->a)
    return(1);
  else
    if (t1->a < t2->a)
      return(-1);
    else
      if (t1->b > t2->b)
	return(1);
      else
	if (t1->b < t2->b)
	  return(-1);
	else
	  if (t1->c > t2->c)
	    return(1);
	  else
	    if (t1->c < t2->c)
	      return(-1);
	    else
	      if (t1->d > t2->d)
		return(1);
	      else
		if (t1->d < t2->d)
		  return(-1);
		else
		  return(0);
}


int compare_angle(angle_rec *t1, angle_rec *t2)
{
  if (t1->a > t2->a)
    return(1);
  else
    if (t1->a < t2->a)
      return(-1);
    else
      if (t1->b > t2->b)
	return(1);
      else
	if (t1->b < t2->b)
	  return(-1);
	else
	  if (t1->c > t2->c)
	    return(1);
	  else
	    if (t1->c < t2->c)
	      return(-1);
		else
		  return(0);
}










