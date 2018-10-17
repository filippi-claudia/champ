/*****
This file is part of the Babel Program
Copyright (C) 1995 Molecular Arts Corporation. All Rights Reserved.

For more information about the M3D file format please contact :
    info@molecules.com
For more information about babel please contact :
    babel@mercury.aichem.arizona.edu
----------------------------------------------------------------------------
FILE : rdm3d.c
AUTHOR(S) : Tony Tribelli <adtribelli@acm.org>
DATE : 12-95
PURPOSE : routines to read a Molecular Arts M3D file
******/

#include "bbltyp.h"

int read_m3d (FILE *file1, ums_type *mol) {
    int     i;
    char    input_line[BUFF_SIZE];
    char    temp_type[8],
            steric[4];
    char    bo_string[16];
    int     column;
  
    fgets (input_line, sizeof (input_line), file1); /* Header */
    fgets (input_line, sizeof (input_line), file1); /* Molecule definition */
    sscanf (input_line, "%d %d",
                        &Atoms,                 /* Number of atoms */
                        &Bonds);                /* Number of bonds */
    ShowProgress (Atoms, "Reading Atoms");
    initialize_ums (&mol);
    steric[0] = '.';                            /* Get steric substring */
    steric[2] = '\0';                           /*  ready to append     */
    column = locate_input_type("M3D");
    for (i = 1; i <= Atoms; i ++) {
        UpdateProgress ();
        fgets (input_line, sizeof (input_line), file1); /* Atom definition */
        sscanf (input_line, "%*d %s %c %lf %lf %lf",
                            temp_type,          /* Element symbol */
                            &steric[1],         /* Steric number */
                            &X (i),             /* X coordinate */
                            &Y (i),             /* Y coordinate */
                            &Z (i));            /* Z coordinate */
        strcat (temp_type, steric);             /* Create element.steric */
    
        Atomic_number(i) = get_input_type (i,column, temp_type, Type (i), dummy);
    }
    for (i = 0; i < Bonds; i++) {
        fgets (input_line, sizeof (input_line), file1); /* Bond definition */
        sscanf (input_line, "%*d %d %d %s",
                            &Start (i),         /* Atom 1 identifier */
                            &End (i),           /* Atom 2 identifier */
                            bo_string);         /* Bond type */
        Bond_order (i) = translate_m3d_bond_order (bo_string);
    }
    dissect_connection_table (mol);
    return (1);
}

int translate_m3d_bond_order (char *bo_string) {
    int ret = 1;                                /* Assume Single */

    switch (*bo_string) {
        case '1' :
            if (*(bo_string+1) == '.')          /* Check for a "1.5" */
                ret = 5;                        /* Resonant/Aromatic */
            break;
        case 'D' :
        case 'd' :
        case '2' :
            ret = 2;                            /* Double */
            break;
        case 'T' :
        case 't' :
        case '3' :
            ret = 3;                            /* Triple */
            break;
        case 'R' :
        case 'r' :
        case 'A' :
        case 'a' :
            ret = 5;                            /* Resonant/Aromatic */
            break;
        case 'H' :
        case 'h' :
        case '0' :
            break;
    }
    return (ret);
}











