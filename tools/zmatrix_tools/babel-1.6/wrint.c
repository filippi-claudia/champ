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

FILE : wrint.c
AUTHOR(S) : Pat Walters
DATE : 11-92
PURPOSE : Routines to write a mopac internal coordinate file, also contains 
internal to cartesian coordinate conversion routines.  This module was adapted
from some of Ajay's original code as well as the routines in MOPAC 5.

******/

#include "bbltyp.h"


int 
  write_mopac_internal(FILE *file1, ums_type *mol)
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

  fprintf(file1,"%s\n",OutputKeywords);
  fprintf(file1,"%s\n",OutfileName);
  fprintf(file1,"%s\n",Title);
  
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
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);    
    fprintf(file1," %-2s%12.7f  %d  %12.6f  %d  %12.6f  %d  %3d  %3d  %3d      %4.4f\n",
	    type_name,
	    mol->internal[i].r,num1,
	    mol->internal[i].w,num2,
	    mol->internal[i].t,num3,
	    mol->internal[i].na,
	    mol->internal[i].nb,
	    mol->internal[i].nc,0.0); 
  }
  return(TRUE);
}


void  
cartint(ums_type *mol)
{
  int i=0; 
  int j=0; 
  int k=0;
  int im1=0;
  double sum=0.0;
  double r=0.0;
  
  mol->internal[2].na=1;

  for(i = 1; i <= Atoms; i++)
  {
    mol->internal[i].na = 2;
    mol->internal[i].nb = 3;
    mol->internal[i].nc = 4;
    im1=i-1; 
    if(im1!=0)
    {
      sum=1.0E30;
      for(j = 1; j <= im1; j++)
      {  
	r = distance(Point(i),Point(j));
	r = SQUARE(r);
	if((r < sum) && (mol->internal[j].na != j) && (mol->internal[j].nb != j))
	{
	  sum = r;
	  k = j;  
	}
      }
      mol->internal[i].na = k;
      
      if(i > 2)
	mol->internal[i].nb = mol->internal[k].na;
      if(i > 3)
	mol->internal[i].nc = mol->internal[k].nb;
    }
  } 
  mol->internal[1].na = 0;
  mol->internal[1].nb = 0;
  mol->internal[1].nc = 0;
  mol->internal[2].nb = 0;
  mol->internal[2].nc = 0;
  mol->internal[3].nc = 0;
}


void 
cartgeom(ums_type *mol)
{  
  int count=0,newcount=0; /* newcount is I1 in the fortran ver */
  int num=0; /* II in the fortran version */
  double TOL = 0.0;
  int j,k,l;
  double angle;
  double sum;
  double r=0.0;


  for(count=2;count <= mol -> num_atoms; count++)
  {
    j=mol->internal[count].na;
    k=mol->internal[count].nb;
    l=mol->internal[count].nc;
    if(count>=3)
    {
      num=count;
      mol->internal[count].w = bond_angle(mol->atoms[num].point,mol->atoms[j].point,mol->atoms[k].point);
      if(count>=4)
      {
	angle=bond_angle(Point(num),Point(j),Point(k));
	       
	TOL = 0.2617994;
	if((angle > (PI - TOL))||(angle < TOL))
	{
	  sum=0.0;
	  for(newcount = 1; newcount <= num-1; newcount++)
	  {
	    r = SQUARE(distance(Point(newcount),Point(k)));

	    if((r < sum)&&(newcount != j)&&(newcount != k))
	    {
	      angle=bond_angle(Point(num),Point(j),Point(k));
	      if((angle < PI-TOL) && (angle > TOL))
	      {
		sum=r;
		l=newcount;
		mol->internal[count].nc = l;		
	      }
	    }
	  }
	  if((sum > 99)&&(TOL > 0.1))
	  {
	    TOL=0.087266;
	    exit(0);
	    
	  }
	}
	mol->internal[count].t = torsion(Point(num),
				    Point(j),
				    Point(k),
				    Point(l));
      }
    }
    mol->internal[count].r = distance(Point(count),Point(j));
  }
  mol->internal[1].r = 0.0;
  mol->internal[1].w = 0.0;
  mol->internal[1].t = 0.0;
  mol->internal[2].w = 0.0;
  mol->internal[2].t = 0.0;
  mol->internal[3].t = 0.0;    
}















