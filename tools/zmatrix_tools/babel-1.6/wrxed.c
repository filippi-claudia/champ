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

FILE : wrxed.c
AUTHOR(S) : Keith Trollope (kit1000@cam.ac.uk)
DATE : 5-95
PURPOSE : Routines to write a XED (new COSMIC) file
******/

#include "bbltyp.h"

int 
write_xed(FILE *file1, ums_type *mol)
{ 
  int i;
  char temp_type[5];
  int result;
  float zero=0.0F;
  int type_name, mass;

  fprintf(file1,"%10.3lf%10i%10i\n",Energy,Atoms,Bonds);
  fprintf(file1,"File conversion by  Babel\n");
  for(i=0;i<Bonds ;i++){
    fprintf(file1,"%8i%8i",Start(i),End(i));
    if ( !((i+1) % 5) ) fputc('\n',file1);
  }
  if (Bonds%5) fputc('\n',file1);
  for(i=1;i<=Atoms;i++) {
    result=xlate_std_type( "XED",Type(i),temp_type);
    if (result == 0)
    {
      fprintf(stderr,"Unable to assign XED type to atom %d type = %s\n",
              i,Type(i));
      type_name = 29;
    }
    else
      type_name = atoi(temp_type);
    switch (type_name) {
           case 1: case 2: case 3: case 4:
             mass=6; break;
           case 5: case 6: case 7: case 8: case 9: case 23: case 25:
             mass=7; break;
           case 10: case 11: case 22: case 24: case 26:
             mass=8; break;
           case 12: case 13:
             mass=16; break;
           case 14:
             mass=15; break;
           case 15:
             mass=1; break;
           case 16:
             mass=9; break;
           case 17:
             mass=17; break;
           case 18:
             mass=35; break;
           case 19:
             mass=53; break;
           default:
             mass=0;
    }
    fprintf(file1,"%6i%15.6lf%15.6lf%15.6lf%6i%12.4f\n",
        mass,X(i),Y(i),Z(i),type_name,zero);
  }
  fprintf(file1,"    1         0.0000    0         0.0000\n");
  return (1);
}

