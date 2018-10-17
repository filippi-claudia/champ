/*****
This file was written as an extension of the Babel Program by:
Mark Benzel
3D Graphics
SGI
mbenzel@sgi.com

The Babel Program info is as follows:
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
------------------------------------------------------------------------------

FILE : wrmiv.c
AUTHOR(S) : Mark Benzel
DATE : 06-96
PURPOSE : Routines to write an Inventor file containing Molecular Inventor
          nodes to be used by Inventor file viewer such as "ivview" from SGI
******/

#include "bbltyp.h"

extern element_type *elements;

static char *level1 = "    ";
static char *level2 = "        ";
static char *level3 = "            ";


int 
write_mol_inventor(FILE *file1, ums_type *mol)
{ 
  int i;
  int includeFile;
  int result;
  int numOnLine;
  char type_name[5];
  char buf[132];
  char *pos;
  FILE *tmpFile;

/*
 * See if the user wants to include all of the original file
 */
  includeFile = 0;
  uppercase(OutputKeywords);
  if (EQ(OutputKeywords, "INCLUDE")) {
    includeFile = 1;
  }

/*
 * Write the Inventor header
 */

  fprintf(file1,"#Inventor V2.1 ascii\n\n");

  fprintf(file1,"# Generated from %s\n\n", InfileName);

  /*
   * Write the root node
   */
  fprintf(file1,"Separator {\n");

  /*
   * Write the ChemData node
   */
  fprintf(file1, "%sChemData {\n", level1);
  fprintf(file1, "%snumberOfAtoms %d\n", level2, Atoms);
  fprintf(file1, "%snumberOfBonds %d\n", level2, Bonds);

  /* Start of the associatedData */
  if ((InfileType != quanta) && (includeFile == 1)) {
    tmpFile = fopen(InfileName, "r");
    if (tmpFile != NULL) {
      fprintf(file1, "%sassociatedData ChemAssociatedData {\n", level2);
      fprintf(file1, "%sdescription \"From %s file %s\"\n", 
        level2, InputTypeName, InfileName);
      fprintf(file1, "%sasciiData [\n", level2);
      while (fgets(buf, 132, tmpFile) != NULL) {
        pos = buf;
        /*
         * Change any " in the input file to ' since " will mess up
         * reading of the file in MolInventor programs.
         */
        while ((pos = strchr(pos, '\"')) != NULL) {
            *pos = '\'';
        }
        strip_return(buf);
        fprintf(file1, "%s\"%s\",\n", level3, buf);
      }
      fclose(tmpFile);
      fprintf(file1, "%s]\n", level3);
      fprintf(file1, "%sbinaryData [ ]\n", level2);
      fprintf(file1, "%s}\n", level2);
    }
  }
  /* End of associatedData */

  /* Start of atomicNumber */
  numOnLine = 0;
  fprintf(file1, "%satomicNumber [\n", level2);
  fprintf(file1, "%s", level3);
  for (i = 1; i <= Atoms; i++) {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1, "%d", get_atomic_number(type_name));
    if (i < Atoms) {
      fprintf(file1, ",");
      if (++numOnLine == 8) {
        fprintf(file1, "\n");
        fprintf(file1, "%s", level3);
        numOnLine = 0;
      }
      else {
        fprintf(file1, " ");
      }
    }
  }
  fprintf(file1, " ]\n");
  /* End of atomicNumber */

  /* Start of atomId */
  numOnLine = 0;
  fprintf(file1, "%satomId [\n", level2);
  fprintf(file1, "%s", level3);
  for (i = 1; i <= Atoms; i++) {
    fprintf(file1, "%d", i);
    if (i < Atoms) {
      fprintf(file1, ",");
      if (++numOnLine == 8) {
        fprintf(file1, "\n");
        fprintf(file1, "%s", level3);
        numOnLine = 0;
      }
      else {
        fprintf(file1, " ");
      }
    }
  }
  fprintf(file1, " ]\n");
  /* End of atomId */

  /* Start of atomName */
  numOnLine = 0;
  fprintf(file1, "%satomName [\n", level2);
  for (i = 1; i <= Atoms; i++) {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    if (i < Atoms) {
      fprintf(file1, "%s\"%s\",\n", level3, type_name);
    }
    else {
      fprintf(file1, "%s\"%s\" ]\n", level3, type_name);
    }
  }
  /* End of atomName */

  /* Start of atomIndex */
  numOnLine = 0;
  fprintf(file1, "%satomIndex [\n", level2);
  fprintf(file1, "%s", level3);
  for (i = 1; i <= Atoms; i++) {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    fprintf(file1, "%d", get_atomic_number(type_name));
    if (i < Atoms) {
      fprintf(file1, ",");
      if (++numOnLine == 8) {
        fprintf(file1, "\n");
        fprintf(file1, "%s", level3);
        numOnLine = 0;
      }
      else {
        fprintf(file1, " ");
      }
    }
  }
  fprintf(file1, " ]\n");
  /* End of atomIndex */

  /* Start of atomCoordinates */
  numOnLine = 0;
  fprintf(file1, "%satomCoordinates [\n", level2);
  for (i = 1; i <= Atoms; i++) {
    if (i < Atoms) {
      fprintf(file1, "%s%f %f %f,\n", level3,
        (float)X(i), (float)Y(i), (float)Z(i));
    }
    else {
      fprintf(file1, "%s%f %f %f ]\n", level3,
        (float)X(i), (float)Y(i), (float)Z(i));
    }
  }
  /* End of atomCoordinates */

  /* Start of bondFrom */
  numOnLine = 0;
  fprintf(file1, "%sbondFrom [\n", level2);
  fprintf(file1, "%s", level3);
  for (i = 0; i < Bonds; i++) {
    fprintf(file1, "%d", Start(i)-1);
    if (i < Bonds - 1) {
      fprintf(file1, ",");
      if (++numOnLine == 8) {
        fprintf(file1, "\n");
        fprintf(file1, "%s", level3);
        numOnLine = 0;
      }
      else {
        fprintf(file1, " ");
      }
    }
  }
  fprintf(file1, " ]\n");
  /* End of bondFrom */

  /* Start of bondTo */
  numOnLine = 0;
  fprintf(file1, "%sbondTo [\n", level2);
  fprintf(file1, "%s", level3);
  for (i = 0; i < Bonds; i++) {
    fprintf(file1, "%d", End(i)-1);
    if (i < Bonds - 1) {
      fprintf(file1, ",");
      if (++numOnLine == 8) {
        fprintf(file1, "\n");
        fprintf(file1, "%s", level3);
        numOnLine = 0;
      }
      else {
        fprintf(file1, " ");
      }
    }
  }
  fprintf(file1, " ]\n");
  /* End of bondTo */

  /* Start of bondType */
  numOnLine = 0;
  fprintf(file1, "%sbondType [\n", level2);
  fprintf(file1, "%s", level3);
  for (i = 0; i < Bonds; i++) {
    if (Bond_order(i) == 3) {
      fprintf(file1, "TRIPLE_BOND");
    }
    else if (Bond_order(i) == 2) {
      fprintf(file1, "DOUBLE_BOND");
    }
    else {
      fprintf(file1, "SINGLE_BOND");
    }
    if (i < Bonds - 1) {
      fprintf(file1, ",");
      if (++numOnLine == 4) {
        fprintf(file1, "\n");
        fprintf(file1, "%s", level3);
        numOnLine = 0;
      }
      else {
        fprintf(file1, " ");
      }
    }
  }
  fprintf(file1, " ]\n");
  /* End of bondType */

  /* Start of bondIndex */
  numOnLine = 0;
  fprintf(file1, "%sbondIndex [\n", level2);
  fprintf(file1, "%s", level3);
  for (i = 0; i < Bonds; i++) {
    fprintf(file1, "%d", i);

    if (i < Bonds - 1) {
      fprintf(file1, ",");
      if (++numOnLine == 8) {
        fprintf(file1, "\n");
        fprintf(file1, "%s", level3);
        numOnLine = 0;
      }
      else {
        fprintf(file1, " ");
      }
    }
  }
  fprintf(file1, " ]\n");
  /* End of bondIndex */

  /*
   * End of ChemData
   */
  fprintf(file1, "%s}\n", level1);

  /*
   * The ChemDisplayParam and ChemUI nodes.  Use the default values.
   */
  fprintf(file1, "%sChemDisplayParam { }\n", level1);
  fprintf(file1, "%sChemUI { }\n", level1);

  /*
   * The ChemRadii node.  Use the radii stored in elements vdw_rad.
   */
  fprintf(file1, "%sChemRadii {\n", level1);
  fprintf(file1, "%satomRadiiBinding RADII_PER_ATOM_INDEXED\n", level2);
  numOnLine = 0;
  fprintf(file1, "%satomRadii [\n", level2);
  fprintf(file1, "%s", level3);
  for (i = 0; i < MAX_ELEMENTS; i++) {
    fprintf(file1, "%f", (float)elements[i].vdw_rad);
    if (i < MAX_ELEMENTS - 1) {
      fprintf(file1, ",");
      if (++numOnLine == 3) {
        fprintf(file1, "\n");
        fprintf(file1, "%s", level3);
        numOnLine = 0;
      }
      else {
        fprintf(file1, " ");
      }
    }
  }
  fprintf(file1, " ]\n");
  fprintf(file1, "%s}\n", level1);

  /*
   * End of ChemRadii
   */

  /*
   * The ChemColor node.  Use the colors stored in elements 
   * red, green, blue.
   */
  fprintf(file1, "%sChemColor {\n", level1);
  numOnLine = 0;
  fprintf(file1, "%satomColorBinding ATOM_PER_ATOM_INDEXED\n", level2);
  fprintf(file1, "%satomColor [\n", level2);
  for (i = 0; i < MAX_ELEMENTS; i++) {
    if (i < MAX_ELEMENTS - 1) {
      fprintf(file1, "%s%f %f %f,\n", level3,
        (float)elements[i].red, (float)elements[i].green,
        (float)elements[i].blue);
    }
    else {
      fprintf(file1, "%s%f %f %f ]\n", level3,
        (float)elements[i].red, (float)elements[i].green,
        (float)elements[i].blue);
    }
  }
  fprintf(file1, "%s}\n", level1);

  /*
   * End of ChemColor
   */

  /*
   * The ChemDisplay.  Use the default values.
   */
  fprintf(file1, "%sChemDisplay { }\n", level1);

  /*
   * End the root node
   */
  fprintf(file1,"}\n");
  return (result);
}
