#include "bbltyp.h"

/*-----------------------------------------------------
Notice of blatant theft.  This code was totally lifted
from Roger Sayle's Internal2Cartesian routine in RasMol 2.6.
Thanks Roger.

Roger's code was a whole lot cleaner than the earlier version
which was stolen from MOPAC 5.  However, by looking at 
Roger's varialbe names it looks like his code was derived
from MOPAC also.  Oh well, what comes around goes around :-)

WPW 4-1-96  April Fools !!!!!
------------------------------------------------------*/

int int_to_cart(ums_type *mol)
{
  register double cosph,sinph,costh,sinth,coskh,sinkh;
  register double cosa,sina,cosd,sind;
  register double dist,angle,dihed;
  
  register double xpd,ypd,zpd,xqd,yqd,zqd;
  register double xa,ya,za,xb,yb,zb;
  register double rbc,xyb,yza,temp;
  register double xpa,ypa,zqa;
  register double xd,yd,zd;
  register int flag;
  int i, na, nb, nc;

  /* Atom #1 */
  X(1) = 0.0;
  Y(1) = 0.0;
  Z(1) = 0.0;
  
  if (Atoms == 1)
    return( TRUE );
  
  /* Atom #2 */
  X(2) = R(2);
  Y(2) = 0.0;
  Z(2) = 0.0;
  
  if( Atoms == 2 )
    return( TRUE );
  
  /* Atom #3 */
  dist = R(3);
  angle = W(3) * DEG_TO_RAD;
  cosa = cos(angle);
  sina = sin(angle);
  
  if( NA(3) == 1 )
  { 
    X(3) = X(1) + cosa*dist;
  } 
  else 
  {   
    X(3) = X(2) - cosa*dist;
  }
  Y(3) = sina*dist;
  Z(3) = 0.0;
  
  for (i = 4; i <= Atoms; i++)
  {   
    dist = R(i);
    angle = W(i) * DEG_TO_RAD;
    dihed = T(i) * DEG_TO_RAD;
    
    na = NA(i);
    nb = NB(i);
    nc = NC(i);
    
    xb = X(nb) - X(na);
    yb = Y(nb) - Y(na);
    zb = Z(nb) - Z(na);
    
    rbc = xb*xb + yb*yb + zb*zb;
    if( rbc < 0.0001 )
      return( FALSE );
    rbc = 1.0/sqrt(rbc);
    
    cosa = cos(angle);
    sina = sin(angle);
    
    
    if( fabs(cosa) >= 0.999999 )
    { 
      /* Colinear */
      temp = dist*rbc*cosa;
      X(i)  = X(na) + temp*xb;
      Y(i) =  Y(na) + temp*yb;
      Z(i) =  Z(na) + temp*zb;
    } 
    else
    {
      xa = X(nc) - X(na);
      ya = Y(nc) - Y(na);
      za = Z(nc) - Z(na);
      
      sind = -sin(dihed);
      cosd = cos(dihed);
      
      xd = dist*cosa;
      yd = dist*sina*cosd;
      zd = dist*sina*sind;
      
      xyb = sqrt(xb*xb + yb*yb);
      if( xyb < 0.1 )
      {  
	/* Rotate about y-axis! */
	temp = za; za = -xa; xa = temp;
	temp = zb; zb = -xb; xb = temp;
	xyb = sqrt(xb*xb + yb*yb);
	flag = TRUE;
      }
      else 
	flag = FALSE;
      
      costh = xb/xyb;
      sinth = yb/xyb;
      xpa = costh*xa + sinth*ya;
      ypa = costh*ya - sinth*xa;
      
      sinph = zb*rbc;
      cosph = sqrt(1.0 - sinph*sinph);
      zqa = cosph*za  - sinph*xpa;
      
      yza = sqrt(ypa*ypa + zqa*zqa);
      
      if( yza > 1.0E-10 )
      {   
	coskh = ypa/yza;
	sinkh = zqa/yza;
	
	ypd = coskh*yd - sinkh*zd;
	zpd = coskh*zd + sinkh*yd;
      } 
      else
      { 
	/* coskh = 1.0; */
	/* sinkh = 0.0; */
	ypd = yd;
	zpd = zd;
      }
      
      xpd = cosph*xd  - sinph*zpd;
      zqd = cosph*zpd + sinph*xd;
      xqd = costh*xpd - sinth*ypd;
      yqd = costh*ypd + sinth*xpd;
      
      if( flag )
      { 
	/* Rotate about y-axis! */
	X(i) = X(na) - zqd;
	Y(i) = Y(na) + yqd;
	Z(i) = Z(na) + xqd;
      } 
      else
      {  
	X(i) = X(na) + xqd;
	Y(i) = Y(na) + yqd;
	Z(i) = Z(na) + zqd;
      }
    }
  }
  return( TRUE );
}





