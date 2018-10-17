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

FILE : rdpsgvbout.c
AUTHOR(S) : Pat Walters, Tom Pollard
DATE : 1-94
PURPOSE : routines to read a PS-GVB output file
******/

#include "bbltyp.h"
#define DELIMS "\n\t "

int
read_psgvb_output(FILE *file1, ums_type *mol)
{
  char the_line[BUFF_SIZE];
  int i;
  int ready;
  long pos = 0;
  long end_pos = 0;
  int result;
  double the_energy = 0;

  while (fgets(the_line,sizeof(the_line), file1) != NULL)
    {
      {
         char *sp;
         for (sp=the_line; *sp; sp++)
             *sp=tolower(*sp);
       }

      if (strstr(the_line,"start of program onee"))
        break;

      if (strstr(the_line,"SCFE"))
        the_energy = atof(gettoken(the_line,DELIMS,5));

      if (   strstr(the_line,"input geometry:") != NULL
          || strstr(the_line,"symmetrized geometry:") != NULL
          || strstr(the_line,"new geometry:") != NULL)
        {
#ifdef DEBUG
          printf("   -> %s\n",the_line);
#endif
          for (i = 0; i < 2; i++)   /* Skip two lines */
            fgets(the_line,sizeof(the_line), file1);

          pos = ftell(file1);       /* Remember where the atoms begin */
          Atoms = 0;
          while (fgets(the_line,sizeof(the_line), file1) != NULL)
            {
              if (count_tokens(the_line,DELIMS) == 0)   /* Blank line */
                break;
              Atoms++;
            }
        }
    }

#ifdef DEBUG
          printf("(done)-> %s\n",the_line);
#endif

  if (feof(file1) && pos == 0)
    {
#ifdef DEBUG
      printf("\nno structure was found\n");
#endif
      Atoms = 0;
      return(FALSE);
    }

  initialize_ums(&mol);
  initialize_internal(&mol);
  Energy = the_energy;

  ShowProgress(Atoms,"Reading Atoms");
  end_pos = ftell(file1);
  fseek(file1,pos,0);
  i = 0;
  while (i < Atoms)
    {
      UpdateProgress();
      fgets(the_line,sizeof(the_line), file1);
      i++;
#ifdef DEBUG
      printf("%d -> %s\n",i,the_line);
#endif
      strcpy(Type(i),gettoken(the_line,DELIMS,1));
      X(i) = atof(gettoken(the_line,DELIMS,2));
      Y(i) = atof(gettoken(the_line,DELIMS,3));
      Z(i) = atof(gettoken(the_line,DELIMS,4));
      clean_atom_type(Type(i));
    }

  if (Atoms > 0)
    {
      result = assign_radii(mol);
      result = assign_bonds(mol);
      result = assign_types(mol);
      result = build_connection_table(mol);
      assign_bond_order(mol);
    }

  fseek(file1,end_pos,0);
/*  read_to_eof(file1);    only read geometry from iteration 1 */
  return(TRUE);

}
