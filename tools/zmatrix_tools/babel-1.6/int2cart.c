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

FILE : int2cart.c
AUTHOR(S) : Pat Walters
DATE : 11-92
PURPOSE : Routines to convert internal to cartesian coordinates
Most of this was adapted from Ajay's routines and those found in MOPAC 5.
******/



#include "bbltyp.h"

#define CONV PI/180.0
static warning wstr;

int old_int_to_cart(ums_type *mol)
{
  
  double *coord[4], *geo[4], *geol[4];
  int *na, *nb, *nc;
  
  int i;
  double ccos,cosa;
  int natoms;
  int ma,mb,mc;
  double xa, ya, za, xb, yb, zb;
  double rbc;
  double xyb;
  int k;
  double yza;
  double xpa, xpb, ypa, xqa, zqa;
  double cosph, sinph, costh, sinth, coskh, sinkh;
  double cosd, sina, sind;
  double xd, yd, zd;
  double xpd, ypd, zpd;
  double xqd, yqd, zqd;
  double xrd;
  

  natoms = Atoms;

    
  for (i = 0; i <= 3; i++)
  {
    coord[i] = (double *)malloc((natoms + 1) * sizeof(double));
    geo[i] = (double *)malloc((natoms + 1) * sizeof(double));
    geol[i] = (double *)malloc((natoms + 1) * sizeof(double));
  }
  
  na = (int *)malloc((natoms + 1) * sizeof(int));
  nb = (int *)malloc((natoms + 1) * sizeof(int));
  nc = (int *)malloc((natoms + 1) * sizeof(int));
  
  for (i = 1; i <= natoms; i++)
  {
    geo[1][i] = mol->internal[i].r;
    geo[2][i] = mol->internal[i].w;
    geo[3][i] = mol->internal[i].t;
    na[i] = mol->internal[i].na;
    nb[i] = mol->internal[i].nb;
    nc[i] = mol->internal[i].nc;
  }
  
  for (i = 1; i <= natoms; i++)
  {
    geol[1][i] = geo[1][i];
    geol[2][i] = geo[2][i] * CONV;
    geol[3][i] = geo[3][i] * CONV;
  }

  if( natoms > 1 ) 
  {
    coord[1][2] = geol[1][2];
    coord[2][2] = 0.0;
    coord[3][2] = 0.0;
  }
  
  if( natoms > 2 ) 
  {
    ccos = cos(geol[2][3]);
    if (na[3] == 1)
      coord[1][3] = coord[1][1] + geol[1][3] * ccos;
    else
      coord[1][3] = coord[1][2] - geol[1][3] * ccos;
    coord[2][3] = geol[1][3] * sin(geol[2][3]);
    coord[3][3] = 0.0;
  }
  
  for (i = 4; i <= natoms; i++)
  {
    cosa = cos(geol[2][i]);
    mb = nb[i];
    mc = na[i];
    xb = coord[1][mb] - coord[1][mc];
    yb = coord[2][mb] - coord[2][mc];
    zb = coord[3][mb] - coord[3][mc];
    rbc = 1.0/sqrt(xb*xb + yb*yb + zb*zb);
    
    if ( fabs(cosa) >= 0.9999999991)
    {
      rbc = geol[1][i] * rbc * cosa;
      coord[1][i] = coord[1][mc] + xb * rbc;
      coord[2][i] = coord[2][mc] + yb * rbc;
      coord[3][i] = coord[3][mc] + zb * rbc;
    }
    else
    {
      ma = nc[i];
      xa = coord[1][ma] - coord[1][mc];
      ya = coord[2][ma] - coord[2][mc];
      za = coord[3][ma] - coord[3][mc];
      
      xyb = sqrt(xb*xb + yb*yb);
      k = -1;
      if (xyb <= 0.10)
      {
	xpa = za;
	za = -xa;
	xa = xpa;
	xpb = zb;
	zb = -xb;
	xb = xpb;
	xyb = sqrt(xb*xb + yb*yb);
	k = 1;
      }
      costh = xb/xyb;
      sinth = yb/xyb;
      xpa = xa * costh + ya * sinth;
      ypa = ya * costh - xa * sinth;
      sinph = zb * rbc;
      cosph = sqrt(fabs(1.00 - sinph * sinph));
      xqa = xpa * cosph + za * sinph;
      zqa = za * cosph - xpa * sinph;
      
      yza = sqrt(ypa*ypa + zqa*zqa);
      if ((yza < 1.0E-1) && (yza > 1.0E-10))
      {
	sprintf(wstr,"THE FAULTY ATOM IS ATOM NUMBER %d",i);
	show_warning(wstr);
	sprintf(wstr,"ATOMS %d, %d AND %d ARE WITHIN %f ANGSTROMS OF A STRAIGHT LINE",mc,mb,ma,yza);
	show_warning(wstr);
      }
      coskh = ypa/yza;
      sinkh = zqa/yza;
      if (yza < 1.0E-10)
      {
	coskh = 1.0;
	sinkh = 0.0;
      }
      
      sina = sin(geol[2][i]);
      sind = -sin(geol[3][i]);
      cosd = cos(geol[3][i]);
      xd = geol[1][i] * cosa;
      yd = geol[1][i] * sina * cosd;
      zd = geol[1][i] * sina * sind;
      
      ypd = yd * coskh - zd * sinkh;
      zpd = zd * coskh + yd * sinkh;
      xpd = xd * cosph - zpd * sinph;
      zqd = zpd * cosph + xd * sinph;
      xqd = xpd * costh - ypd * sinth;
      yqd = ypd * costh + xpd * sinth;
      if (k >= 1)
      {
	xrd = -zqd;
	zqd = xqd;
	xqd = xrd;
      }
      coord[1][i] = xqd + coord[1][mc];
      coord[2][i] = yqd + coord[2][mc];
      coord[3][i] = zqd + coord[3][mc];
    }
  }
  for (i = 1; i <= natoms; i++)
  {
    X(i) = coord[1][i];
    Y(i) = coord[2][i];
    Z(i) = coord[3][i];
  }
  
  for (i = 0; i <= 3; i++)
  {
    free(coord[i]);
    free(geo[i]);
    free(geol[i]);
  }
  free(na);
  free(nb);
  free(nc);
  return(TRUE);
}

