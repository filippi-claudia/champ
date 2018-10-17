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

FILE : torlist.c
AUTHOR(S) : Pat Walters
DATE : 3-20-96
PURPOSE : list heavy atom torsions

******/

#include "bbltyp.h"

#define IsHeavyAtm(x)  (Atomic_number(x) > 1)

int print_torsion_list(FILE *file1,ums_type *mol) 
{
  int a,b,c,d;
  int i,j,k;
  double angle1;
  int angle_count = 0;
  torsion_rec *tr;

  tr = (torsion_rec *)malloc(Atoms * 10 * sizeof(torsion_rec));
  if (tr == NULL)
  {
    printf("Memory Allocation Error at %s %s \n",__FILE__,__LINE__);
    return(FALSE);
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
	    if (IsHeavyAtm(a) && IsHeavyAtm(b) && IsHeavyAtm(c) && IsHeavyAtm(d))
	    {
	      tr[angle_count].a = a;
	      tr[angle_count].b = b;
	      tr[angle_count].c = c;
	      tr[angle_count].d = d;
	      angle_count ++;
	    }
	  }
      }
  }

  qsort(tr,angle_count,sizeof(torsion_rec),QSORT_PROTO compare_torsion);
  for (i = 0; i < angle_count; i++)
  {
    angle1 = torsion(Point(tr[i].a),Point(tr[i].b),Point(tr[i].c),Point(tr[i].d));
    fprintf(file1,"%10.3f ",angle1);
  }

  fprintf(file1,"%10.3f\n",Energy);
  free(tr); 
}	




