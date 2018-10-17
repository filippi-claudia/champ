/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : rdpdb.c
AUTHOR(S) : Pat Walters
DATE : 10-92
PURPOSE : Routines to read a Brookhave Protien Data Bank file

******/
#include "bbltyp.h"

#define DELIMS "\t\n "
#define PDBATOM 0
#define PDBHETATM 1

#include    <unistd.h>
#include    <sys/time.h>             
#include    <sys/times.h>
#include    <sys/resource.h>

int 
read_pdb(FILE *file1, ums_type *mol)
{
  char pdb_line[BUFF_SIZE];
  int i = 0;
  int result;
  char temp_char;
  long int pos;
  int chain_num = 1;
  char res_num_str[5];
  char the_atm_id[5],the_res_name[5]; 
  char the_serial_num[6];

#ifdef DOTIME
  struct timeval s_val;
  double startTime, endTime;
#endif

  pos = ftell(file1);

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  startTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
#endif

  while ((fgets(pdb_line,sizeof(pdb_line), file1) != NULL) && NOTEQn(pdb_line,"END",3))
  { 
    if (EQn(pdb_line,"ATOM",4) || EQn(pdb_line,"HETATM",6))
    {
      i++;
    }
  }

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  endTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
  printf("Getting Atoms: %lf\n", endTime - startTime);
  startTime = endTime;
#endif

  fseek(file1,pos,0);
  Atoms = i;

  initialize_ums(&mol);
  initialize_residues(&mol);
  ShowProgress(Atoms,"Reading Atoms");
	 
  i = MIN_ATOM;
  while ((fgets(pdb_line,sizeof(pdb_line), file1) != NULL) && NOTEQn(pdb_line,"END",3))
  {
    if (EQn(pdb_line,"TER",3))
      chain_num++;
    if (EQn(pdb_line,"ATOM",4) || EQn(pdb_line,"HETATM",6))
    {
      UpdateProgress();
      parse_atom_record(&pdb_line[6], i, mol);
      if (pdb_line[0] == 'A')
        stringToType(PDBATOM, i, mol);
      else 
        stringToType(PDBHETATM, i, mol);
      ChainNum(i) = chain_num;
      i ++;
    }
  }

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  endTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
  printf("Reading Atoms: %lf\n", endTime - startTime);
  startTime = endTime;
#endif

  /*
   * Make sure Valence is zeroed out
   */

  Bonds = 0;
  for (i = MIN_ATOM; i <= Atoms; i++)
  {
    Valence(i) = 0;
  }

  fseek(file1,pos,0);
  process_connect_records(file1,mol); 

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  endTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
  printf("Reading CONECT: %lf\n", endTime - startTime);
  startTime = endTime;
#endif

  result = assign_radii(mol);  
  if (!result) 
  {
    release_ums(mol);
    return(FALSE);
  }

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  endTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
  printf("Assigning Radii: %lf\n", endTime - startTime);
  startTime = endTime;
#endif

  if ((mol->control) &&  (NOTEQ(InputKeywords,"CONECT_ONLY")))
  {
    result = assign_pdb_bonds(mol);
    if (!result)
      {
	release_ums(mol);
	return(FALSE);
      }
  }
  

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  endTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
  printf("Assigning PDB bonds: %lf\n", endTime - startTime);
  startTime = endTime;
#endif

  qsort(mol->connections,Bonds,sizeof(connect_type),QSORT_PROTO sort_connections);
  dissect_connection_table(mol);

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  endTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
  printf("Dissect Connection Table: %lf\n", endTime - startTime);
  startTime = endTime;
#endif

  if ((mol->control) &&  (NOTEQ(InputKeywords,"CONECT_ONLY")))
    check_bonds(mol);

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  endTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
  printf("Check Bonds: %lf\n", endTime - startTime);
  startTime = endTime;
#endif

  result = assign_types(mol);

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  endTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
  printf("Assigning Types: %lf\n", endTime - startTime);
  startTime = endTime;
#endif

  if (!result)
  {
    release_ums(mol);
    return(FALSE);
  }
  assign_bond_order(mol);

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  endTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
  printf("Assigning Bond Order: %lf\n", endTime - startTime);
  startTime = endTime;
#endif

  result = build_connection_table(mol);

#ifdef DOTIME
  gettimeofday(&s_val, 0);
  endTime = (double) s_val.tv_sec + 0.000001*s_val.tv_usec;
  printf("Building Connection Table: %lf\n", endTime - startTime);
  startTime = endTime;
#endif

  if (!result)
  {
    release_ums(mol);
    return(FALSE);
  }
  qsort(mol->connections,Bonds,sizeof(connect_type),QSORT_PROTO sort_connections);
  dissect_connection_table(mol);
  BO(0,0) = -1;
  return(TRUE);  
}

void fix_A_type(char *type, char *id, char *res_type)
{
  /* Check for ASP, ASN and ASX */
  /* Check for GLU, GLN and GLX */
  if (EQn(id,"AD1",3) && (EQn(res_type,"AS",2) || EQn(res_type,"N",1)))
    strcpy(type,"O");
  else
    if (EQn(id,"AD2",3) && (EQn(res_type,"AS",2) || EQn(res_type,"N",1)))
      strcpy(type,"N");
  else
    if (EQn(id,"AE1",3) && (EQn(res_type,"GL",2) || EQn(res_type,"Q",1)))
      strcpy(type,"O");
  else
    if (EQn(id,"AE2",3) && (EQn(res_type,"GL",2) || EQn(res_type,"Q",1)))
      strcpy(type,"N");
  else
    if (EQn(id,"AD1",3) && (EQn(res_type,"HIS",3) || EQn(res_type,"H",1)))
      strcpy(type,"N");
  else
    if (EQn(id,"AE1",3) && (EQn(res_type,"HIS",3) || EQn(res_type,"H",1)))
      strcpy(type,"C");
  else
    if (EQn(id,"AE2",3) && (EQn(res_type,"HIS",3) || EQn(res_type,"H",1)))
      strcpy(type,"N");
  else
    if (EQn(id,"AD2",3) && (EQn(res_type,"HIS",3) || EQn(res_type,"H",1)))
      strcpy(type,"C");
}

void process_connect_records(FILE *file1, ums_type *mol)
{
  char pdb_line[BUFF_SIZE];  
  int the_serial_num, conn1, conn2, conn3, conn4;

  while ((fgets(pdb_line,sizeof(pdb_line), file1) != NULL) && NOTEQn(pdb_line,"END",3))
  {
    if EQn(pdb_line,"CONECT",6)
    {
      parse_conect_record(&pdb_line[6], &the_serial_num, &conn1, &conn2,
        &conn3, &conn4);
      if (conn1 > 0) add_bond(mol, the_serial_num, conn1);
      if (conn2 > 0) add_bond(mol, the_serial_num, conn2);
      if (conn3 > 0) add_bond(mol, the_serial_num, conn3);
      if (conn4 > 0) add_bond(mol, the_serial_num, conn4);
    }
  }
}

void add_bond(ums_type *mol, int i, int j)
{
  int k, start, end;
  /*
   * Need to find the atoms with the passed in serial numbers
   */

  if (i < j)
  {
    start = i;
    end = j;
  }
  else
  {
    start = j;
    end = i;
  }

  for (k = 1; k <= Atoms; k++)
  {
    if (SerialNum(k) == start)
    {
      start = k;
      break;
    }
  }
  for (k = start + 1; k <= Atoms; k++)
  {
    if (SerialNum(k) == end)
    {
      end = k;
      break;
    }
  }
  if (ChainNum(start) == ChainNum(end))
  {
    Connection(start,Valence(start)) = end;
    Connection(end,Valence(end)) = start;
    Valence(start)++;
    Valence(end)++;
    Start(Bonds) = start;
    End(Bonds) = end;
    Bonds++;
  }
}

int parse_atom_record(char *pdb_line, int the_atom, ums_type *mol)
{
/* ATOMFORMAT "(i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,1x,i3)" */

  int len;
  int sc;
  char tmp_str[10];
  char *tmp_ptr;

  len = strlen(pdb_line) - 1;
  sc = 0;

  /* serial number */
  my_strncpy(tmp_str, &pdb_line[sc], 5); sc += 6;
  SerialNum(the_atom) = atoi(tmp_str);
  if (sc >= len) return(FALSE);

  /* atom name */
  my_strncpy(tmp_str, &pdb_line[sc], 4); sc += 4;
  strcpy(AtmId(the_atom), tmp_str);
  tmp_ptr = rtrim(tmp_str);
  strcpy(Type(the_atom), tmp_str);
  if (sc >= len) return(FALSE);

  /* alternate location  NOTUSED */
  my_strncpy(tmp_str, &pdb_line[sc], 1); sc += 1;
  if (sc >= len) return(FALSE);

  /* residue name */
  my_strncpy(tmp_str, &pdb_line[sc], 3); sc += 4;
  if ((tmp_str[0] == ' ') && (tmp_str[1] == ' ') && (tmp_str[2] == ' '))
  {
    strcpy(ResName(the_atom), "UNK");
  }
  else
  {
    strcpy(ResName(the_atom), tmp_str);
  }
  if (sc >= len) return(FALSE);

  /* chain ID  NOTUSED */
  my_strncpy(tmp_str, &pdb_line[sc], 1); sc += 1;
  if (sc >= len) return(FALSE);

  /* residue sequence number */
  my_strncpy(tmp_str, &pdb_line[sc], 4); sc += 4;
  ResNum(the_atom) = atoi(tmp_str);
  if (sc >= len) return(FALSE);

  /* residue insertion code NOTUSED */
  my_strncpy(tmp_str, &pdb_line[sc], 1); sc += 4;
  if (sc >= len) return(FALSE);

  /* X, Y, Z */
  my_strncpy(tmp_str, &pdb_line[sc], 8); sc += 8;
  X(the_atom) = atof(tmp_str);
  if (sc >= len) return(FALSE);

  my_strncpy(tmp_str, &pdb_line[sc], 8); sc += 8;
  Y(the_atom) = atof(tmp_str);
  if (sc >= len) return(FALSE);

  my_strncpy(tmp_str, &pdb_line[sc], 8); sc += 8;
  Z(the_atom) = atof(tmp_str);
  if (sc >= len) return(TRUE);

#ifdef COMPLETE_PDB
  /* occupancy NOTUSED */
  my_strncpy(tmp_str, &pdb_line[sc], 6); sc += 6;
  if (sc >= len) return(TRUE);

  /* temperature factor NOTUSED */
  my_strncpy(tmp_str, &pdb_line[sc], 6); sc += 7;
  if (sc >= len) return(TRUE);

  /* footnote number NOTUSED*/
  my_strncpy(tmp_str, &pdb_line[sc], 3);
#endif

  return(TRUE);
}

int
parse_conect_record(char *pdb_line, int *serial_number,
  int *conn1, int *conn2, int *conn3, int *conn4)
{
/* CONECTFORMAT "(5i5)" */
  int len;
  int sc;
  char tmp_str[10];

  len = strlen(pdb_line) - 1;
  sc = 0;
             
  *serial_number = *conn1 = *conn2 = *conn3 = *conn4 = 0;

  my_strncpy(tmp_str, &pdb_line[sc], 5); sc += 5;
  *serial_number = atoi(tmp_str);
  if (sc >= len) return TRUE;
      
  my_strncpy(tmp_str, &pdb_line[sc], 5); sc += 5;
  *conn1 = atoi(tmp_str);
  if (sc >= len) return TRUE;
         
  my_strncpy(tmp_str, &pdb_line[sc], 5); sc += 5;
  *conn2 = atoi(tmp_str);
  if (sc >= len) return TRUE;
    
  my_strncpy(tmp_str, &pdb_line[sc], 5); sc += 5;
  *conn3 = atoi(tmp_str);
  if (sc >= len) return TRUE;
      
  my_strncpy(tmp_str, &pdb_line[sc], 5); sc += 5;
  *conn4 = atoi(tmp_str);
      
  return TRUE;  
}

int
stringToType(int record_type, int the_atom, ums_type *mol)
{
   int ilen;
   int atnum;
   char *tmpstr1;
   char *tmpstr2;
   char *pos;
         
   if ((ilen = (int)strlen ( Type(the_atom) )) == 0) {
      show_warning ( "Error parsing the atom name (length of zero)." );
      return (FALSE);
   }
            
   tmpstr1 = strdup ( Type(the_atom) );
   tmpstr2 = strdup ( Type(the_atom) );
   uppercase ( tmpstr1 );
               
   if (record_type == PDBATOM) {
      /*   
       * In the standard PDB files, the first character will be blank or
       * a digit.  The second character will be the element we are looking
       * for.  There are exceptions to this.  Often D is used for deuterium
       * so we will convert D to H.  For some of the residues, the atom type
       * cannot be determined and the element name will start with A.
       * We will look for these and handle as appropriate.
       */
      ljust ( tmpstr2 );
      ljust ( tmpstr1 );
      pos = tmpstr1;
      while (isdigit(pos[0])) pos++;
      pos[1] = '\0';
      if (pos[0] == 'A') {
         fix_A_type(Type(the_atom), tmpstr2, ResName(the_atom));
         free ( tmpstr1 );
         free ( tmpstr2 );
         return (TRUE);
      }
   } else {  
      /*
       * The HETATM record is more complicated than the ATOM record in terms
       * of determining the type of element from the atom label.
       * We will do the following:
       * If there is no residue name then we assume the atom label is the
       * name of the element.   
       * If the residue name and the atom label are the same then we will
       * assume that the atom label is the name of the element.
       * Generally, if there is a letter in the first column then the
       * atom label represents an element with a 2 letter element name.
       * The exceptions so far to this is when the residue name is
       * ADR, COA, FAD, GPG, NAD, NAL or NDP.  In these cases there is a
       * letter in the first column and the second letter is the element
       * name.
       */  
      pos = tmpstr1;
      /*
       * Check the residue name and if one of those listed above, move
       * the pointer.
       */
      if ((!strcmp ( ResName(the_atom), "ADR" )) ||
          (!strcmp ( ResName(the_atom), "COA" )) ||
          (!strcmp ( ResName(the_atom), "FAD" )) ||
          (!strcmp ( ResName(the_atom), "GPG" )) ||
          (!strcmp ( ResName(the_atom), "NAD" )) ||
          (!strcmp ( ResName(the_atom), "NAL" )) ||
          (!strcmp ( ResName(the_atom), "NDP" ))) {
         pos++;
         pos[1] = '\0';
      } else if (isdigit(pos[0])) {
         pos++;
         pos[1] = '\0';   
      } else if (pos[0] != ' ') {
         pos[2] = '\0';
         if (isalpha(pos[1]) && isupper(pos[1])) pos[1] = tolower(pos[1]);
      } else {
         ljust ( pos );
         pos[1] = '\0';
      }
   }   
   strcpy(Type(the_atom), pos);
   free ( tmpstr1 );
   free ( tmpstr2 );
   return (TRUE);
}

#undef DELIMS
#undef PDBATOM
#undef PDBHETATM
