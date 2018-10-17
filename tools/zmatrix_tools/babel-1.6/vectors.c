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

FILE : vectors.c 
AUTHOR(S) : Matt Stahl
DATE : 6-5-93
PURPOSE : Various vector manipulation routines

Modified 6-20-93 by Pat Walters
Added matrix routines

******/

#include "bbltyp.h"

void point_to_vect(coord_type pt, vect_type *v)
{
  v->x = pt.x;
  v->y = pt.y;
  v->z = pt.z;
}

void pts_2_vect(ums_type *mol, vect_type *vect,int pt1,int pt2)
{
  
  vect->x = X(pt1) - X(pt2);
  vect->y = Y(pt1) - Y(pt2);
  vect->z = Z(pt1) - Z(pt2);
}


double vect_ang(vect_type *vect1,vect_type *vect2)
{
  double angle;
  double mag;
  double dp;
  double rad_ang;

  mag = magnitude(vect1) * magnitude(vect2);
  dp = dot(vect1,vect2)/mag;
  if (dp < -0.999999)
    dp = -0.9999999;
  if (dp > 0.9999999)
    dp = 0.9999999;
  rad_ang = acos(dp);
  angle = RAD_TO_DEG * rad_ang;
  return(angle);
}


double dot(vect_type *vect1,vect_type *vect2)
{
  double dot_prod;
  
  dot_prod = vect1->x*vect2->x + vect1->y*vect2->y + vect1->z*vect2->z;
  return(dot_prod);
}

void cross_prod(vect_type *vect1,vect_type *vect2,vect_type *normal)
{
  normal->x = vect1->y*vect2->z - vect2->y*vect1->z;
  normal->y = vect2->x*vect1->z - vect1->x*vect2->z;
  normal->z = vect1->x*vect2->y - vect2->x*vect1->y;
}

double magnitude(vect_type *vect1)
{
  double mag;
  
  mag = sqrt(SQUARE(vect1->x) + SQUARE(vect1->y) + SQUARE(vect1->z));

  return(mag);
}

void normalize_vect(vect_type *v1)
{
  double mag;
  
  mag = magnitude(v1);

  if (mag != 0.0)
  {
    v1->x = v1->x/mag;
    v1->y = v1->y/mag;
    v1->z = v1->z/mag;
  }
}

void vect_sum(vect_type *vect1, vect_type *vect2, vect_type *vect_sm)
{
  vect_sm->x = vect1->x + vect2->x;
  vect_sm->y = vect1->y + vect2->y;
  vect_sm->z = vect1->z + vect2->z;
}

void vect_diff(vect_type *vect1, vect_type *vect2, vect_type *vect_sm)
{
  vect_sm->x = vect1->x - vect2->x;
  vect_sm->y = vect1->y - vect2->y;
  vect_sm->z = vect1->z - vect2->z;
}

void scal_x_vect(vect_type *vect1, float scalar)
{
   vect1->x = (vect1->x * scalar);
   vect1->y = (vect1->y * scalar);
   vect1->z = (vect1->z * scalar);
}

coord_type point_plus_vector(coord_type *p1, vect_type *v1)
{
  coord_type result;
  
  result.x = p1->x + v1->x;
  result.y = p1->y + v1->y;
  result.z = p1->z + v1->z;

  return(result);
}

coord_type point_times_vector(coord_type *p1, vect_type *v1)
{
  coord_type result;
  
  result.x = p1->x * v1->x;
  result.y = p1->y * v1->y;
  result.z = p1->z * v1->z;

  return(result);
}

double determinant_3x3(matrix_3x3 *m)
{

  double x,y,z;

  x = m->a1 * (m->b2 * m->c3 - m->b3 * m->c2);
  y = m->a2 * (m->b3 * m->c1 - m->b1 * m->c3);
  z = m->a3 * (m->b1 * m->c2 - m->b2 * m->c1);

  return(x + y + z);
}

void invert_vector(vect_type *v)
{
  v->x *= -1.0;
  v->y *= -1.0;
  v->z *= -1.0;
}

void invert_3x3(matrix_3x3 *m)
{
  matrix_3x3 t;
  double det;

  det = determinant_3x3(m);

  if (det != 0.0)
  {
    t.a1 = m->b2*m->c3 - m->b3*m->c2;
    t.b1 = m->b3*m->c1 - m->b1*m->c3;
    t.c1 = m->b1*m->c2 - m->b2*m->c1;
    t.a2 = m->c2*m->a3 - m->c3*m->a2;
    t.b2 = m->c3*m->a1 - m->c1*m->a3;
    t.c2 = m->c1*m->a2 - m->c2*m->a1;
    t.a3 = m->a2*m->b3 - m->a3*m->b2;
    t.b3 = m->a3*m->b1 - m->a1*m->b3;
    t.c3 = m->a1*m->b2 - m->a2*m->b1;
    
    m->a1 = t.a1/det;
    m->b1 = t.b1/det;
    m->c1 = t.c1/det;
    m->a2 = t.a2/det;
    m->b2 = t.b2/det;
    m->c2 = t.c2/det;
    m->a3 = t.a3/det;
    m->b3 = t.b3/det;
    m->c3 = t.c3/det;
  }
}

void dump_3x3(matrix_3x3 *m)
{
  printf("%10.4f %10.4f %10.4f \n",m->a1,m->b1,m->c1);
  printf("%10.4f %10.4f %10.4f \n",m->a2,m->b2,m->c2);
  printf("%10.4f %10.4f %10.4f \n",m->a3,m->b3,m->c3);
}


void mat_3x3_dot_vect(matrix_3x3 *m, vect_type *v)
{
  vect_type temp;
  
  temp.x = m->a1 * v->x + m->b1 * v->y + m->c1 * v->z;
  temp.y = m->a2 * v->x + m->b2 * v->y + m->c2 * v->z;
  temp.z = m->a3 * v->x + m->b3 * v->y + m->c3 * v->z;

  v->x = temp.x;
  v->y = temp.y;
  v->z = temp.z;
}







