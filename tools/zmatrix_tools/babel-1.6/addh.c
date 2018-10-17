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
FILE : addh.c
AUTHOR(S) : Pat Walters
DATE : 5-18-93
PURPOSE : Routines to add hydrogen to a structure
Techniques based on
M. Nardelli
Computers in Chemistry
Vol 6., No. 3, pp 139-152, 1982

Modified 6-18-93 to include Matt's
routines from vectors.h
******/

#include "bbltyp.h"

static warning wstr;

void add_hydrogens(ums_type *mol)
{
  int old_count, num_H_to_add;
  int result;
  
  old_count = Atoms;

  if ((Bond_order(1) == 0) || (BO(0,0) == -1))
    { 
      num_H_to_add = count_missing_hydrogens(mol);   
      place_hydrogens1(mol,old_count,num_H_to_add);
    }
  else
    {
      num_H_to_add = count_missing_bo_hydrogens(mol); 
      place_hydrogens2(mol,old_count,num_H_to_add);
    }
  
#if 0
  if ((mol->control) && (Verbose))
    {
      sprintf(wstr,"the structure is missing %d hydrogens",to_add);
      show_warning(wstr);
    }
#endif
}

void place_hydrogens1(ums_type *mol, int old_count, int num_H_to_add)
{
  int h_count, result;
  int i;
  char temp_title[BUFF_SIZE];

  Atoms += num_H_to_add;
  h_count = old_count + 1;
  strcpy(temp_title,Title);
  result = reinitialize_ums(&mol);
  zero_out_ums(mol,old_count + 1);
  strcpy(Title,temp_title);
  
  for (i = 1; i <= old_count; i++)
    {
      if EQ(Type(i),"C3")
	     switch(Valence(i))
	       {
	       case 3 :
		 add_tertiary_hydrogen(mol,i,h_count,SP3_C_H_DIST);
		 h_count += 1;
		 break;
	       case 2 :
		 add_methylene_hydrogens(mol,i,h_count,SP3_C_H_DIST);
		 h_count += 2;
		 break;
	       case 1 :
		 add_methyl_hydrogen(mol,i,h_count,SP3_C_H_DIST);
		 h_count += 1;
		 add_methylene_hydrogens(mol,i,h_count,SP3_C_H_DIST);
		 h_count += 2;
		 break;
	       }
      else
	if EQ(Type(i),"N3+")
	       switch(Valence(i))
		 {
		 case 2 :
		   add_methylene_hydrogens(mol,i,h_count,SP3_N_H_DIST);
		   h_count += 2;
		   break;
		 case 1 :
		   add_methyl_hydrogen(mol,i,h_count,SP3_N_H_DIST);
		   h_count += 1;
		   add_methylene_hydrogens(mol,i,h_count,SP3_N_H_DIST);
		   h_count += 2;
		   break;
		 }
	else
	  if (EQ(Type(i),"C2") || EQ(Type(i),"Car"))
	    switch(Valence(i))
	      {
	      case 2 :
		add_sp2_hydrogen(mol,i,h_count,SP2_C_H_DIST);
		h_count += 1;
		break;
	      }
	  else
	    if (EQ(Type(i),"Npl") || EQ(Type(i),"Nam") || EQ(Type(i),"Ng+"))
	      switch(Valence(i))
		{
		case 2 :
		  add_sp2_hydrogen(mol,i,h_count,SP2_N_H_DIST);
		  h_count += 1;
		  break;
		}
	    else
	      if EQ(Type(i),"C1")
	{
	  if (Valence(i) == 1)
	    {
	      add_sp_hydrogen(mol,i,h_count,SP_C_H_DIST);
	      h_count++;
	    }
	}
	      else
		if EQ(Type(i),"O3")
	{
	  if (Valence(i) == 1)
	    {
	      add_methyl_hydrogen(mol,i,h_count,SP3_O_H_DIST);
	      h_count++;
	    }
	}
    }
  
  for (i = 1; i <= old_count; i++)
    {
      if EQ(Type(i),"C2")
	     switch(Valence(i))
	       {
	       case 1 :
		 add_vinyl_hydrogens(mol,i,h_count,SP2_C_H_DIST);
		 h_count += 2;
		 break;
	       }
      if (EQ(Type(i),"Npl") || EQ(Type(i),"Nam") || EQ(Type(i),"Ng+"))
	switch(Valence(i))
	  {
	  case 1 :
	    add_vinyl_hydrogens(mol,i,h_count,SP2_C_H_DIST);
	    h_count += 2;
	    break;
	  }
    }
  Atoms = h_count - 1;
  build_connection_table(mol);
}

void add_methyl_hydrogen(ums_type *mol, int c_num, int h_num, double b_length)
{
  int c1,c2,c3;
  vect_type v;
  
  c1 = c_num;
  c2 = Connection(c1,0);
  if (Connection(c2,0) != c1)
    c3 = Connection(c2,0);
  else
    c3 = Connection(c2,1);

  pts_2_vect(mol,&v,c2,c3);
  normalize_vect(&v);
  scal_x_vect(&v,(float) b_length);
  Point(h_num) = point_plus_vector(&Point(c1),&v);

  type_added_hydrogen(mol,c1,h_num);
}


void add_sp_hydrogen(ums_type *mol, int c_num, int h_num, double b_length)
{
  int c1, c2;
  vect_type v;
  
  c1 = c_num;
  c2 = Connection(c1,0);
  pts_2_vect(mol,&v,c1,c2);
  normalize_vect(&v);
  scal_x_vect(&v,(float) b_length);
  Point(h_num) = point_plus_vector(&Point(c1),&v);
  
  type_added_hydrogen(mol,c1,h_num);
}

void add_sp2_hydrogen(ums_type *mol, int c_num, int h_num, double b_length)
{
  int c1, c2, c3;  
  vect_type v,v1,s;
  

  c1 = c_num;
  c2 = Connection(c1,0);
  c3 = Connection(c1,1);

  pts_2_vect(mol,&v,c1,c2);
  normalize_vect(&v);
  pts_2_vect(mol,&v1,c1,c3);
  normalize_vect(&v1);
  vect_sum(&v,&v1,&s);
  normalize_vect(&s);
  scal_x_vect(&s,(float) b_length);
  Point(h_num) = point_plus_vector(&Point(c1),&s);

  type_added_hydrogen(mol,c1,h_num);
}

void add_vinyl_hydrogens(ums_type *mol, int c_num, int h_num, double b_length)
{
  int c1,c2,c3,c4;
  int h1,h2;
  vect_type v,v1;
  
  h1 = h_num;
  h2 = h1 + 1;

  c1 = c_num;
  c2 = Connection(c1,0);
  if (Connection(c2,0) != c1)
    c3 = Connection(c2,0);
  else
    c3 = Connection(c2,1);
  if ((Connection(c2,0) != c1) && (Connection(c2,0) != c3))
    c4 = Connection(c2,0);
  else
    if ((Connection(c2,1) != c1) && (Connection(c2,1) != c3))
      c4 = Connection(c2,1);
    else
      c4 = Connection(c2,2);

  pts_2_vect(mol,&v,c2,c3);
  normalize_vect(&v);
  pts_2_vect(mol,&v1,c2,c4);
  normalize_vect(&v1);
  scal_x_vect(&v,(float) b_length);
  Point(h1) = point_plus_vector(&Point(c1),&v);
  scal_x_vect(&v1,(float) b_length);
  Point(h2) = point_plus_vector(&Point(c1),&v1);
  
  type_added_hydrogen(mol,c1,h1);
  type_added_hydrogen(mol,c1,h2);
}

void add_sp3_N_hydrogen(ums_type *mol, int n_num, int h_num, double b_length)
{
  int c1, c2, c3;  
  int h1;
  vect_type v,v1,n,s,h1vect;

  c1 = n_num;
  c2 = Connection(c1,0);
  c3 = Connection(c1,1);
  h1 = h_num;

  pts_2_vect(mol,&v,c1,c2);
  normalize_vect(&v);
  pts_2_vect(mol,&v1,c1,c3);
  normalize_vect(&v1);
  vect_sum(&v,&v1,&s);
  cross_prod(&v,&v1,&n);
  scal_x_vect(&s,(float) ONE_OVER_SQRT3);
  scal_x_vect(&n,(float) SQRT_TWO_THIRDS);
  vect_sum(&s,&n,&h1vect);
  normalize_vect(&h1vect);
  scal_x_vect(&h1vect,(float) b_length);
  Point(h1) = point_plus_vector(&Point(c1),&h1vect);

  type_added_hydrogen(mol,c1,h1);
}
  
void add_methylene_hydrogens(ums_type *mol, int c_num, int h_num, double b_length)
{
  int c1, c2, c3;  
  int h1, h2;
  vect_type v,v1,n,s,h1vect,h2vect;

  c1 = c_num;
  c2 = Connection(c1,0);
  c3 = Connection(c1,1);
  h1 = h_num;
  h2 = h_num+1;

  pts_2_vect(mol,&v,c1,c2);
  normalize_vect(&v);
  pts_2_vect(mol,&v1,c1,c3);
  normalize_vect(&v1);
  vect_sum(&v,&v1,&s);
  cross_prod(&v,&v1,&n);
  scal_x_vect(&s,(float) ONE_OVER_SQRT3);
  scal_x_vect(&n,(float) SQRT_TWO_THIRDS);
  vect_sum(&s,&n,&h1vect);
  normalize_vect(&h1vect);
  scal_x_vect(&h1vect,(float) b_length);
  Point(h1) = point_plus_vector(&Point(c1),&h1vect);

  scal_x_vect(&n,-1.0);
  vect_sum(&s,&n,&h2vect);
  normalize_vect(&h2vect);
  scal_x_vect(&h2vect,(float) b_length);
  Point(h2) = point_plus_vector(&Point(c1),&h2vect);
  
  type_added_hydrogen(mol,c1,h1);
  type_added_hydrogen(mol,c1,h2);
}

void add_tertiary_hydrogen(ums_type *mol, int c_num, int h_num, double b_length)
{
  vect_type v1,v2,v3,s;
  int c1,c2,c3,c4,h1;
  matrix_3x3 m;

  c1 = c_num;
  c2 = Connection(c1,0);
  c3 = Connection(c1,1);
  c4 = Connection(c1,2);
  h1 = h_num;
  
  
  pts_2_vect(mol,&v1,c1,c2);
  normalize_vect(&v1);
  pts_2_vect(mol,&v2,c1,c3);
  normalize_vect(&v2);
  pts_2_vect(mol,&v3,c1,c4);
  normalize_vect(&v3);
  
  m.a1 = v1.x;  m.b1 = v1.y; m.c1 = v1.z;
  m.a2 = v2.x;  m.b2 = v2.y; m.c2 = v2.z;
  m.a3 = v3.x;  m.b3 = v3.y; m.c3 = v3.z;
  
  invert_3x3(&m);

  s.x = m.a1 + m.b1 + m.c1;
  s.y = m.a2 + m.b2 + m.c2;
  s.z = m.a3 + m.b3 + m.c3;

  normalize_vect(&s);
  scal_x_vect(&s,(float) b_length);  
  Point(h1) = point_plus_vector(&Point(c1),&s);
  type_added_hydrogen(mol,c1,h1);
}

int
type_added_hydrogen(ums_type *mol,int c1, int h1)
{
  
  Atomic_number(h1) = 1;
  if (Atomic_number(c1) == 6)
    strcpy(Type(h1),"HC");
  else
    strcpy(Type(h1),"H");

  if (HasResidues)
  {
    strcpy(AtmId(h1),"H");
    strcpy(ResName(h1),ResName(c1));
    ResNum(h1) = ResNum(c1);
  }
  
  Valence(h1) = 1;
  Connection(c1,Valence(c1)) = h1;
  BO(c1,Valence(c1)) = 1;
  Valence(c1)++;
  Connection(h1,0) = c1;
  BO(h1,0) = 1;
  return(TRUE);
}


int
count_missing_hydrogens(ums_type *mol)
{
  int missing = 0;
  int i;
  char temp_type[5];
  int type_valence;
  int result;

  for (i = 1; i <= Atoms; i++)
    {
      result = xlate_std_type("HAD",Type(i),temp_type);
      if (result == 0)
	{
	  sprintf(wstr,"Unable to assign valence to atom %d type = %s",
		  i,Type(i));
	  show_warning(wstr);
	  strcpy(temp_type,"0");
	}
      type_valence = atoi(temp_type);
      if ((Valence(i) < type_valence) && (Valence(i) > 0))
	{
	  missing += type_valence - Valence(i);
	}
#if 0
      printf("num = %d type = %s val = %d type_valence = %d\n",
	     i,Type(i),Valence(i),type_valence); 
#endif
    }
  return(missing);
}


void
add_2d_hydrogens(ums_type *mol)
{

  int old_count, to_add, h_count;
  int type_valence;
  char temp_type[10];
  int result;
  int i,j;

  old_count = Atoms;
  to_add = count_missing_hydrogens(mol); 

  Atoms += to_add;
  h_count = old_count + 1;
  result = reinitialize_ums(&mol);
  zero_out_ums(mol,old_count + 1);

  for (i = 1; i <= old_count; i++)
    {
      to_add = 0;

      result = xlate_std_type("HAD",Type(i),temp_type);
      if (result == 0)
	{
	  sprintf(wstr,"Unable to assign valence to atom %d type = %s",
		  i,Type(i));
	  show_warning(wstr);
	  strcpy(temp_type,"0");
	}
      type_valence = atoi(temp_type);
      if ((Valence(i) < type_valence) && (Valence(i) > 0))
	{
	  to_add = type_valence - Valence(i);
	  /*       printf("num = %d type = %s val = %d \n",i,Type(i),Valence(i)); */
	}

    
      for(j = 0;j < to_add;j++)
	{
	  Connection(i,Valence(i)) = h_count;
	  Valence(i)++;
	  Start(Bonds) = i;
	  End(Bonds) = h_count;
	  Bond_order(Bonds) = 1;
	  Bonds++;
	  sprintf(Type(h_count),"H%c",Type(i)[0]);
	  Connection(h_count,0) = i;
	  Valence(h_count)++;
	  h_count++;
	}
    }
  
}












