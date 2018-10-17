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

FILE : gastchg.c
AUTHOR(S) : Pat Walters
DATE : 10-95
PURPOSE : calculate Gasteiger-Marsilli charges
******/

#include "bbltyp.h"

#define NUM_SIGMA_PARAM 14

void calc_gasteiger_charges(ums_type *mol)
{
  calc_sigma_charges(mol);
}

void setup_sigma_params(gast_param *par)
{
  int i;
  
  strcpy(par[0].atm_typ,  "H"); par[0].a =   7.17; par[0].b =   6.24; par[0].c =  -0.56;
  strcpy(par[1].atm_typ, "C3"); par[1].a =   7.98; par[1].b =   9.18; par[1].c =   1.88;
  strcpy(par[2].atm_typ, "C2"); par[2].a =   8.79; par[2].b =   9.32; par[2].c =   1.51;
  strcpy(par[3].atm_typ, "C1"); par[3].a =  10.39; par[3].b =   9.45; par[3].c =   0.73;
  strcpy(par[4].atm_typ, "N3"); par[4].a =  11.54; par[4].b =  10.82; par[4].c =   1.36;
  strcpy(par[5].atm_typ, "N2"); par[5].a =  12.87; par[5].b =  11.15; par[5].c =   0.85;
  strcpy(par[6].atm_typ, "N1"); par[6].a =  15.68; par[6].b =  11.70; par[6].c =  -0.27;
  strcpy(par[7].atm_typ, "O3"); par[7].a =  14.18; par[7].b =  12.92; par[7].c =   1.39;
  strcpy(par[8].atm_typ, "O2"); par[8].a =  17.07; par[8].b =  13.79; par[8].c =   0.47;
  strcpy(par[9].atm_typ, "F");  par[9].a =  14.66; par[9].b =  13.85; par[9].c =   2.31;
  strcpy(par[10].atm_typ,"Cl"); par[10].a = 11.00; par[10].b =  9.69; par[10].c =  1.35;
  strcpy(par[11].atm_typ,"Br"); par[11].a = 10.08; par[11].b =  8.47; par[11].c =  1.16;
  strcpy(par[12].atm_typ,"I");  par[12].a =  9.90; par[12].b =  7.96; par[12].c =  0.96;
  strcpy(par[13].atm_typ,"S3");  par[13].a = 10.14; par[13].b =  9.13; par[13].c =  1.38;

  for (i = 0; i < NUM_SIGMA_PARAM; i++)
    par[i].d = par[i].a + par[i].b + par[i].c;
}
                              
void lookup_sigma_params(ums_type *mol, gast_param *master_par, gast_param *atm_par)
{
  int i,j;
  int found;
  char type_name[5];
  
  for (i = 1; i <= Atoms; i++)
  {
    get_output_type(i,"MAP",Type(i),type_name,all_caps);
    for (j = 0; j < NUM_SIGMA_PARAM; j++)
    {
      found = FALSE;
      if EQn(type_name,master_par[j].atm_typ,2)
      {
        found = TRUE;
        strcpy(atm_par[i].atm_typ,master_par[j].atm_typ);
        atm_par[i].a = master_par[j].a;
        atm_par[i].b = master_par[j].b;
        atm_par[i].c = master_par[j].c;
        atm_par[i].d = master_par[j].d;
        break;
      }
    }
    if (!found)
    {
      strcpy(atm_par[i].atm_typ,"Du");
      atm_par[i].a = 0.0;
      atm_par[i].b = 0.0;
      atm_par[i].c = 0.0;
      atm_par[i].d = 1.0;
    }
  }
}


void calc_sigma_charges(ums_type *mol)
{
  gast_param mast_par[NUM_SIGMA_PARAM], *atm_par;
  int bonds;
  double *xx;
  double d1,d2,z1,q1;
  int i,j,t;
  int cycle = 0;
  
#ifdef DEBUG
  input_chg = (double *)malloc((Atoms + 1) * sizeof(double));
  for (i = 1; i <= Atoms; i++)
    input_chg[i] = Charge(i);
#endif
  
  atm_par = (gast_param *)malloc((Atoms + 1) * sizeof(gast_param));
  xx = (double *)malloc((Atoms + 1) * sizeof(double));
  
  setup_sigma_params(mast_par);
  lookup_sigma_params(mol,mast_par,atm_par);

  /* set formal charges */

  for (i = 1; i <= Atoms; i++)
  {
    switch(Atomic_number(i))
    {
    case 6 :
      if (EQ(Type(i),"C+"))
	Charge(i) = 1.0;
      break;
    case 7 :
      bonds = count_attached_bonds(mol,i);
      if (bonds == 4) 
      {
	Charge(i) = 1.0;
      }
	
      if (EQ(Type(i),"N3+"))
	Charge(i) = 1.0;
      break;
    case 8 :
      if (EQ(Type(i),"O-"))
	Charge(i) = -0.5;
      break;
    default :
      Charge(i) = 0.0;
    }

    xx[i] = atm_par[i].a;
  }
  z1 = 1.0;
  
  do
  {
    z1 *= 0.5;
    d1 = 0.0;
    
    for (i = 1; i <= Atoms; i++)
      if (atm_par[i].d != 1.0)
    {
      q1 = Charge(i);
      for (j = 0; j < Valence(i); j++)
        if (atm_par[Connection(i,j)].d != 1.0)
      {
        t = Connection(i,j);
        d2 = atm_par[t].d;
        if (xx[t] > xx[i])
          d2 = atm_par[i].d;
        
        if (Type(i)[0] == 'H')
          d2 = 20.02;
        if (Type(t)[0] == 'H')
          d2 = 20.02;
        
        Charge(i) = Charge(i) + (xx[t] - xx[i])/d2*z1;


      }
      q1 = fabs(Charge(i) - q1);
      if (q1 > d1)
        d1 = q1;
    }
    
    if (d1 >= 0.001)
      for (i = 1; i <= Atoms; i++)
      {
        xx[i] = atm_par[i].a + atm_par[i].b * Charge(i) + atm_par[i].c * SQUARE(Charge(i));
      }
    
    cycle++;
    
  } while ((d1  > 0.001) && (cycle <= 5));

#ifdef DEBUG
  for (i = 1; i <= Atoms; i++)
    printf("%2d %4s %10.4f %10.4f %10.4f\n",
           i,Type(i),Charge(i),input_chg[i], fabs(Charge(i) - input_chg[i]));
  free(input_chg);
#endif

  free(atm_par);
  free(xx);
}

void print_gasteiger_params(ums_type *mol, gast_param *atm_par)
{
  int i;
  
  for (i = 1; i <= Atoms; i++)
  {
    printf("%4s %10.3f %10.3f %10.3f %10.3f\n",Type(i),atm_par[i].a,atm_par[i].b,atm_par[i].c,atm_par[i].d);
  }
}

int is_carboxylate(ums_type *mol, int atm)
{
  int i,conn;
  
  for (i = 0; i < Valence(atm); i++)
  {
    conn = Connection(atm,i);
    if (EQ(Type(conn),"Cac"))
      return(TRUE);
  }
    return(FALSE);
}















