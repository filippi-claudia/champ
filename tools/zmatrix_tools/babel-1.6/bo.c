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

FILE : bndord.c
AUTHOR(S) : Pat Walters
DATE : 2-10-93
PURPOSE : Assign bond orders based on atom types, clean up conjugated pi systems
******/

#include "bbltyp.h"

static char wstr[256];
static int cycles;

#define SINGLE_DOUBLE_CUTOFF 0.95
#define DOUBLE_TRIPLE_CUTOFF 0.81


void assign_bond_order2(ums_type *mol)
{
  int *dots;
  int i,j;
  bnd_stack stk;
  set_type *dbatoms, *dbbonds;
  ring_struct rings;
  int conn;
  
  for (i = 1; i <= Atoms; i++)
    Redo(i) = 0;

  dots = (int *)malloc((Atoms + 1) * sizeof(int));
  stk.bond = (int *)malloc((Bonds) * sizeof(int));
  stk.choice = (int *)malloc((Bonds) * sizeof(int));

  dbatoms = init_set_minbits(Atoms);
  dbbonds = init_set_minbits(Bonds);

  find_SSSR(mol,&rings);
  tag_ring_atoms(mol,rings.ring_list,rings.count);

  assign_hybrid_radii(mol);
  estimate_bond_order2(mol);

  for (i = 0; i < rings.count; i++)
  {
    if (rings.ring_list[i].length == 5)
      process_5_ring(mol,&rings,i);
  }  

  dissect_connection_table(mol);

  cleanup_rings(&rings);

  
  for (i = 0; i < Bonds; i++)
  {
    if (Bond_order(i) == 2)
    {
      if ((Valence(Start(i)) > 1) && (Type(Start(i))[0] == 'O'))
	Bond_order(i) = 1;
      else
	if ((Valence(End(i)) > 1) && (Type(End(i))[0] == 'O'))
	  Bond_order(i) = 1;
    }
  }

  for (i = 0; i < Bonds; i++)
  {
    if (Bond_order(i) == 2)
    {
      if ((Valence(Start(i)) == 1) && (Type(Start(i))[0] == 'N'))
	Bond_order(i) = 1;
      else
	if ((Valence(End(i)) == 1) && (Type(End(i))[0] == 'N'))
	  Bond_order(i) = 1;
    }
  }

  for (i = 1; i <= Atoms; i++)
    dots[i] = 0;

  for (i = 0; i < Bonds; i++)
    if (Bond_order(i) > 1)
    {
      if ((Redo(Start(i))) && (check_for_carbonyl(mol,Start(i))) != 3)
	dots[Start(i)] = 1;
      if ((Redo(End(i))) &&(check_for_carbonyl(mol,End(i))) != 3)
	dots[End(i)] = 1;

      if ((Valence(Start(i)) == 1) || (Valence(End(i)) == 1))
      {
	dots[Start(i)] = 0;
	dots[End(i)] = 0;
      }
    }

  for (i = 1; i <= Atoms; i++)
  {
    if (EQ(Type(i),"Npl") && (Valence(i) == 3))
      dots[i] = 0;
  }

  stk.ptr = 0;
  cycles = 0;

  connect_the_dots(mol,1,0,dots,&stk);

  for (i = 1; i <= stk.ptr; i++)
  {
    if (Bond_order(stk.bond[i]) > 1)
    {
      biton(dbbonds,stk.bond[i]);
      biton(dbatoms,Start(stk.bond[i]));
      biton(dbatoms,End(stk.bond[i]));
    }
  }

  for (i = 0; i < Bonds; i++)
  {
    if (EQ(Type(Start(i)),"O2") || EQ(Type(End(i)),"O2"))
    {
      biton(dbbonds,i);
      biton(dbatoms,Start(i));
      biton(dbatoms,End(i));
    }
    else
      if (EQ(Type(Start(i)),"O-") && (Valence(Start(i)) == 1))
      {
      biton(dbbonds,i);
      biton(dbatoms,Start(i));
      biton(dbatoms,End(i));
      }
      else
	if (EQ(Type(End(i)),"O-") && (Valence(End(i)) == 1))
	{
	  biton(dbbonds,i);
	  biton(dbatoms,Start(i));
	  biton(dbatoms,End(i));
	}
  }

#ifdef DEBUG
  setprint(dbatoms,"dbatoms");
  setprint(dbbonds,"dbbonds");
#endif

  for (i = 1; i <= Atoms; i++)
    dots[i] = 0;

  for (i = 0; i < Bonds; i++)
  {
    if (Bond_order(i) > 1) 
    {

      if (!bit_is_on(dbatoms,Start(i)) && !bit_is_on(dbatoms,End(i))) 
      {
	dots[Start(i)] = 1;
	dots[End(i)] = 1;
      }
    }
  }

  stk.ptr = 0;
  connect_the_dots(mol,1,0,dots,&stk);

  for (i = 1; i <= stk.ptr; i++)
  {
      biton(dbbonds,stk.bond[i]);
  }


  for (i = 0; i < Bonds; i++)
  {
    if (!bit_is_on(dbbonds,i))
      Bond_order(i) = 1;
  }

  dissect_connection_table(mol);
  
  for (i = 0; i <= Atoms; i++)
  {
    Redo(i) = 0;
    for (j = 0; j < Valence(i); j++)
      Redo(i) += BO(i,j);

    if ((Atomic_number(i) == 6 || Atomic_number(i) == 7) && (Redo(i) > 4))
    {
      for (j = 0; j < Valence(i); j++)
      {
	conn = -1;
	if (BO(i,j) == 2)
	{
	  BO(i,j) = 1;
	  conn = Connection(i,j);
	  break;
	}
      }
      
      if (conn > 0)
      {
	for (j = 0; j < Valence(conn); j++)
	{
	  if (Connection(conn,j) == i)
	    BO(conn,j) = 1;
	}
      }
    }
  }

  build_connection_table(mol);
  
  if (dots)
    free(dots);
  if (stk.bond)
    free(stk.bond);
  if (stk.choice)
    free(stk.choice);
  free_set(dbatoms);
  free_set(dbbonds);
}


void set_bo(ums_type *mol)
{
  int i;
  int sum_code;

  for (i = 0; i < Bonds; i++)
  {
    sum_code = assign_bond_code(Type(Start(i))) +
      assign_bond_code(Type(End(i)));
    switch(sum_code)
    {
    case 6 :
      Bond_order(i) = 3;
      break;
    case 4 :
      Bond_order(i) = 2;
      break;
    default :
      Bond_order(i) = 1;
    }
    if (is_carboxyl(mol,i))
      Bond_order(i) = 2;
    if (Bond_order(i) < 1 || Bond_order(i) > 3)
    {
      sprintf(wstr,"Bond %d - atoms %d - %d is wierd - Bond order is %d\n",
	     i,Start(i),End(i),Bond_order(i));
    show_warning(wstr);
    }
  }
}

int connect_the_dots(ums_type *mol, int atm, int start, int *dots, bnd_stack *stk)
{
  int i;
  int con;
  int bond;
  int done;
  int new_atm, choice_bnd;

  if (start == Valence(atm))
    return(0);

  if (dots[atm])
  {
    done = FALSE;
    for (i = start; (i < Valence(atm) && !done); i++)
    {
      con = Connection(atm,i);
      if (dots[con])
      {
	stk->ptr++;
	bond = get_bond_number(mol,atm,con);
	stk->bond[stk->ptr]= bond;
	if (atm == Start(bond))
	  stk->choice[stk->ptr] = i+1;
	else
	  stk->choice[stk->ptr] = -(i+1);
	dots[atm] -= 1;
	dots[con] -= 1;
	done = TRUE;
      }
    }
    
    if ((!done) && (stk->ptr > 0))
    {
      bond = stk->bond[stk->ptr];
      choice_bnd = stk->choice[stk->ptr];
      if (choice_bnd > 0)
	    new_atm = Start(bond);
      else
	new_atm = End(bond);
      choice_bnd = abs(choice_bnd);
      stk->ptr--;
      dots[Start(bond)] += 1;
      dots[End(bond)] += 1;
      connect_the_dots(mol,new_atm,choice_bnd,dots,stk);
    }
  }

  if (cycles > 10000)
    return(0);

  if (atm == Atoms) 
    return(1);
  else
  {
    cycles++;
    connect_the_dots(mol,atm+1,0,dots,stk);
  }
  return(TRUE);
}

int get_bond_number(ums_type *mol,int start, int end)
{
  int i;
  
  for (i = 0; i < Bonds; i++)
  {
    if ((Start(i) == start) && (End(i) == end))
      return(i);
    else
      if ((Start(i) == end) && (End(i) == start))
	return(i);
  }
  return(-1);
}


void tag_ring_atoms(ums_type *mol, path *ring_set, int ring_count)
{
  int i,j;

  for (i = 1; i <= Atoms; i++)
    Redo(i) = 0;

  for (i = 0; i < ring_count; i++)
  {
    if (ring_set[i].bogus == FALSE)
    {
      for (j = 0; j < ring_set[i].length; j++)
      {
	Redo(ring_set[i].path_atoms[j]) = 1;
      }
    }
  }
}

double get_bond_ratio(ums_type *mol, int a1, int a2)
{
  double dist,cov_sum;
  
  dist = distance(Point(a1),Point(a2));
  cov_sum = BORadius(a1) + BORadius(a2);
  return(dist/cov_sum);
}


void estimate_bond_order2(ums_type *mol)
{
  double ratio;
  int bo;
  int i;
  char start_type[5];
  char end_type[5];
  
  for (i = 0; i < Bonds; i++)
  {
    bo = 1;
    ratio = get_bond_ratio(mol,Start(i),End(i));
    get_output_type(Start(i),"HYB",Type(Start(i)),start_type,all_caps);
    get_output_type(End(i),"HYB",Type(End(i)),end_type,all_caps);

    if (ratio <= DOUBLE_TRIPLE_CUTOFF)
    {
      if ((start_type[0] == '1') && (end_type[0] == '1'))
	bo = 3;
    }
    else
      if (ratio <= SINGLE_DOUBLE_CUTOFF)
      {
	get_output_type(Start(i),"HYB",Type(Start(i)),start_type,all_caps);
	get_output_type(End(i),"HYB",Type(End(i)),end_type,all_caps);
	if ((start_type[0] == '2') && (end_type[0] == '2'))
	  bo = 2;
      }
#ifdef DEBUG
    printf("i = %d Start = %d End = %d r1 = %10.3f r2 = %10.3f dist = %10.2f sum = %10.2f ratio = %10.2f bond order = %4d\n",
	   i,Start(i),End(i),BORadius(Start(i)),BORadius(End(i)),dist,cov_sum,ratio,bo);
#endif
    Bond_order(i) = bo;
  }
}



void process_5_ring(ums_type *mol, ring_struct *rings, int num)
{
  int bond;
  int a1,a2,a3,a4,a5;
  double t1,t2,t3,t4,t5;
  double ratio;

  a1 = rings->ring_list[num].path_atoms[0];
  a2 = rings->ring_list[num].path_atoms[1];
  a3 = rings->ring_list[num].path_atoms[2];
  a4 = rings->ring_list[num].path_atoms[3];
  a5 = rings->ring_list[num].path_atoms[4];

  t1 = torsion(Point(a5),Point(a1),Point(a2),Point(a3));
  t2 = torsion(Point(a1),Point(a2),Point(a3),Point(a4));
  t3 = torsion(Point(a2),Point(a3),Point(a4),Point(a5));
  t4 = torsion(Point(a3),Point(a4),Point(a5),Point(a1));
  t5 = torsion(Point(a4),Point(a5),Point(a1),Point(a2));

#if 0
  printf("%4d %4d %4d %4d %4d\n",a1,a2,a3,a4,a5);
  printf("%10.3f %10.3f %10.3f %10.3f %10.3f\n",t1,t2,t3,t4,t5);
  printf("%10.3f %10.3f %10.3f %10.3f %10.3f\n",T1,T2,T3,T4,T5);
#endif

  if (fabs(t1) < 7.0)
  {
    bond = get_bond_number(mol,a1,a2);
    Bond_order(bond) = 1;
    ratio = get_bond_ratio(mol,a1,a2);
    if (ratio < SINGLE_DOUBLE_CUTOFF)
      Bond_order(bond) = 2;
  }
  if (fabs(t2) < 7.0)
  {
    bond = get_bond_number(mol,a2,a3);
    Bond_order(bond) = 1;
    ratio = get_bond_ratio(mol,a2,a3);
    if (ratio < SINGLE_DOUBLE_CUTOFF)
      Bond_order(bond) = 2;
  }
  if (fabs(t3) < 7.0)
  {
    bond = get_bond_number(mol,a3,a4);
    Bond_order(bond) = 1;
    ratio = get_bond_ratio(mol,a3,a4);
    if (ratio < SINGLE_DOUBLE_CUTOFF)
      Bond_order(bond) = 2;
  }
  if (fabs(t3) < 7.0)
  {
    bond = get_bond_number(mol,a4,a5);
    Bond_order(bond) = 1;
    ratio = get_bond_ratio(mol,a4,a5);
    if (ratio < SINGLE_DOUBLE_CUTOFF)
      Bond_order(bond) = 2;
  }
  if (fabs(t3) < 7.0)
  {
    bond = get_bond_number(mol,a5,a1);
    Bond_order(bond) = 1;
    ratio = get_bond_ratio(mol,a5,a1);
    if (ratio < SINGLE_DOUBLE_CUTOFF)
      Bond_order(bond) = 2;
  }
}



/*-----------------------------------------------
Testing functions for validating bond order assignments
------------------------------------------------*/

void print_bo(ums_type *mol)
{
  int i;
  
  for (i = 0; i < Bonds; i++)
    printf("%4d - %4d %4d %4d\n",i,Start(i),End(i),Bond_order(i));
}

void reset_bonds(ums_type *mol, int *temp_bo)
{
  int i;
  
  for (i = 0; i < Bonds; i++)
    Bond_order(i) = temp_bo[i];
}


void save_bond_orders(ums_type *mol, int *temp_bo)
{
  int i;
  
  for (i = 0; i < Bonds; i++)
    temp_bo[i] = Bond_order(i);
}


int check_bond_order(ums_type *mol)
{
  int mdl_count[4] = {0,0,0,0};
  int babel_count[4] = {0,0,0,0};
  int i;

  for (i = 0; i < Bonds; i++)
  {
    mdl_count[Bond_order(i)]++;
    Bond_order(i) = 1;
  }

  assign_bond_order2(mol);

  for (i = 0; i < Bonds; i++)
  {
    babel_count[Bond_order(i)]++;
  }

  for (i = 1; i <= 3; i++)
  {
    if (mdl_count[i] != babel_count[i])
      return(FALSE);
  }
  return(TRUE);
}



double torsion_angle(coord_type a, coord_type b, coord_type c, coord_type d)
{
  vect_type c1,c2,c3,c4;
  vect_type v1,v2,v3, p,q;
  double xtheta, theta, absth, s;
  int done = FALSE;
  
  point_to_vect(a,&c1);
  point_to_vect(b,&c2);
  point_to_vect(c,&c3);
  point_to_vect(d,&c4);

  vect_diff(&c1,&c2,&v1);
  vect_diff(&c2,&c3,&v2);
  vect_diff(&c3,&c4,&v3);
  
  cross_prod(&v2,&v1,&p);
  cross_prod(&v3,&v2,&q);

  normalize_vect(&p);
  normalize_vect(&q);

  xtheta = dot(&p,&q);
  
  if (xtheta > 1.0)
    xtheta = 1.0;
  if (xtheta < -1.0)
    xtheta = -1.0;
  theta = acos(xtheta);
  
  theta *= 57.29578;
  
  absth = fabs(theta);
  if (absth < 0.001)
  {
    done = TRUE;
    theta = 0.0;
  }
  else if (fabs(absth - 180.0) < 0.001)
  {
    done = TRUE;
    theta = 180.0;
  }
  
  if (!done)
  {
    s = dot(&v1,&q);
    if (s < 0.0)
      theta = 360.0 - theta;
  }
  
  if (theta > 180.0)
    theta -= 360.0;

  return(theta);
}

/*-------------------------------------
dearomatize - turn a ums with aromatic bonds
into a ums with single and double bonds

Inputs:

mol - a ums
-------------------------------------*/

void dearomatize(ums_type *mol)
{
  int *dots;
  bnd_stack stk;  
  int i;

  dots = (int *)malloc((Atoms + 1) * sizeof(int));
  stk.bond = (int *)malloc((Bonds) * sizeof(int));
  stk.choice = (int *)malloc((Bonds) * sizeof(int));

  for (i = 1; i <= Atoms; i++)
    dots[i] = 0;
  
  for (i = 0; i < Bonds; i++)
  {
    if (Bond_order(i) == 5)
    {
      Bond_order(i) = 1;
      dots[Start(i)] = 1;
      dots[End(i)] = 1;
    }
  }

/*  
  for (i = 1; i <= Atoms; i++)
  {
    if (EQ(Type(i),"Car") || EQ(Type(i),"Nar"))
      dots[i] = 1;
    else
      dots[i] = 0;
  }

  for (i = 1; i <= Atoms; i++)
    printf("%d %s %d\n",i,Type(i),dots[i]);
*/
  stk.ptr = 0;
  connect_the_dots(mol,1,0,dots,&stk);

  for (i = 1; i <= stk.ptr; i++)
  {
    Bond_order(stk.bond[i]) = 2;
  }

  dissect_connection_table(mol);

  if (dots)
    free(dots);
  if (stk.bond)
    free(stk.bond);
  if (stk.choice)
    free(stk.choice);
}


int has_aromatic_bonds(ums_type *mol)
{
  int i;
  
  for (i = 0; i < Bonds; i++)
  {
    if (Bond_order(i) == 5)
      return(TRUE);
  }
  return(FALSE);
}

      


























