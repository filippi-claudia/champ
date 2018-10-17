/*****
This file is part of the Babel Program
Copyright (C) 1995 Molecular Arts Corporation. All Rights Reserved.

For more information about the M3D file format please contact :
    info@molecules.com
For more information about babel please contact :
    babel@mercury.aichem.arizona.edu
------------------------------------------------------------------------------
FILE : wrm3d.c
AUTHOR(S) : Tony Tribelli <adtribelli@acm.org>
DATE : 12-95
PURPOSE : Write a Molecular Arts M3D file
******/

#include "bbltyp.h"

int write_m3d (FILE *file1, ums_type *mol) {
    int     i;
    char    temp_type[8],
            *steric;
    int     result;
    char    bond_string[16];
  
    fprintf (file1, "%-16.16s %2d.%2.2d %5d\n",
                    "STRUCTURE",                /* File type */
                    1,                          /* Major version */
                    2,                          /* Minor version */
                    1);                         /* Number of structures */
    fprintf (file1, "%5d %5d % 8.3f %5d %.50s\n",
                    Atoms,                      /* Number of atoms */
                    Bonds,                      /* Number of bonds */
                    0.0F,                       /* Molecular charge */
                    0,                          /* Molecule status */
                    "");                        /* Molecular label */
    for (i = 1; i <= Atoms; i++) {
        result = get_output_type (i, "M3D", Type (i), temp_type, dummy);
        steric = strchr (temp_type, '.');       /* Should have a seperator */
        if (steric != NULL)                     /* Terminate element string */
            *(steric++) = '\0';                 /*  and reference steric    */
        else
            steric = "0";                       /* Default to "none" */
        fprintf (file1, "%5d %-2.2s %c % 8.3f % 8.3f % 8.3f % 6.3f %.30s\n",
                        i,                      /* Atom identifier */
                        temp_type,              /* Element symbol */
                        *steric,                /* Steric number */
                        X (i),                  /* X coordinate */
                        Y (i),                  /* Y coordinate */
                        Z (i),                  /* Z coordinate */
                        0.0F,                   /* Atom partial charge */
                        "");                    /* Atom label */
    }
    for (i = 0; i < Bonds; i++) {
        switch (Bond_order (i)) {
            case 2 :
                strcpy (bond_string, "Double");
                break;
            case 3 :
                strcpy (bond_string, "Triple");
                break;
            case 5 :
                strcpy (bond_string, "Resonant");
                break;
            default :
                strcpy (bond_string, "Single");
        }
        fprintf (file1, "%5d %5d %5d %-.20s\n",
                        i + 1,                  /* Bond identifier */
                        Start (i),              /* Atom 1 identifier */
                        End (i),                /* Atom 2 identifier */
                        bond_string);           /* Bond type */
    }
    return (1);
}
