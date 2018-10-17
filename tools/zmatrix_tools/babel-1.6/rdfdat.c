/*****
  This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
  
  For more information please contact :
  
  babel@mercury.aichem.arizona.edu
  ------------------------------------------------------------------------
  FILE : rdfdat.c
  AUTHOR(S) : Pat Walters
  DATE : 12-93
  PURPOSE : routines to read a CSD FDAT file
******/

#include "bbltyp.h"
#undef PRINT_HEADER
#undef DEBUG

static warning wstr;

int read_fdat(FILE *file1, ums_type *mol)
{
  char csd_line[BUFF_SIZE];
  fract_type f;
  int remark_lines;
  int remark_chars;
  int symmetry_lines;
  long int pos;
  
  char REFCODE[9];
  char NAT[4], NSAT[4], NCON[4], CELL[2];
  char NRFAC[4],NREM[3],NDIS[4],NERR[4],NOPR[4],NRAD[4];
  char A[7],B[7],C[7],ALPHA[7],BETA[7],GAMMA[7];
  char Ap[2],Bp[2],Cp[2],ALPHAp[2],BETAp[2],GAMMAp[2];
  int i,j,k;
  int result;
  matrix_3x3 m;
  int pending = FALSE;
  int ncon;
  
  char previous[15];
  int tokens, start;
  char Xstr[10], Ystr[10], Zstr[10];
  int done = FALSE;
  
  strcpy(previous,"XXXXXXXXXX");
  uppercase(InputKeywords);

  while ((fgets(csd_line,sizeof(csd_line), file1) != NULL) && (!done))
  {
    if (check_for_eof(file1))
      return(FALSE);
    if (csd_line[0] == '#')
    {
      remark_lines = 0;
      done = TRUE;
      my_strncpy(REFCODE,&csd_line[1],8);
      if (check_refcode(InputKeywords,csd_line,mol)) 
      {
	pending = TRUE;
	strcpy(REFCODE,gettoken(REFCODE," ",1));
	my_strncpy(NRFAC,&csd_line[26],3);
	my_strncpy(NREM,&csd_line[29],3);
	my_strncpy(NDIS,&csd_line[32],3);
	my_strncpy(NERR,&csd_line[35],3);
	my_strncpy(NOPR,&csd_line[38],3);
	my_strncpy(NRAD,&csd_line[41],3);
	my_strncpy(NAT,&csd_line[44],3);
	my_strncpy(NSAT,&csd_line[47],3);
	my_strncpy(NCON,&csd_line[53],3);
	my_strncpy(CELL,&csd_line[56],3);
	Atoms = atoi(NAT) + atoi(NSAT);
	ShowProgress(Atoms,"Reading Atoms");
	result = initialize_ums(&mol);
	strcpy(Title,REFCODE);
#ifdef PRINT_HEADER
	printf("REFCODE = %s  NAT = %s NSAT = %s Atoms = %d\n",
	       REFCODE,NAT,NSAT,Atoms);
	printf("NRFAC = %s NREM = %s NDIS = %s NERR =%s NOPR = %s NRAD = %s\n",
	       NRFAC,NREM,NDIS,NERR,NOPR,NRAD);
#endif
	if ((CELL[0] == '1') && (Atoms > 0))     /* Make sure we have cell params */
	{
	  fgets(csd_line,sizeof(csd_line), file1);
	  my_strncpy(A,&csd_line[0],6);
	  my_strncpy(B,&csd_line[6],6);
	  my_strncpy(C,&csd_line[12],6);
	  my_strncpy(ALPHA,&csd_line[18],6);
	  my_strncpy(BETA,&csd_line[24],6);
	  my_strncpy(GAMMA,&csd_line[30],6);
	  my_strncpy(Ap,&csd_line[36],1);
	  my_strncpy(Bp,&csd_line[37],1);
	  my_strncpy(Cp,&csd_line[38],1);
	  my_strncpy(ALPHAp,&csd_line[39],1);
	  my_strncpy(BETAp,&csd_line[40],1);
	  my_strncpy(GAMMAp,&csd_line[41],1);
	  f.A = my_atof(A)/pow(10.0,my_atof(Ap));
	  f.B = my_atof(B)/pow(10.0,my_atof(Bp));
	  f.C = my_atof(C)/pow(10.0,my_atof(Cp));
	  f.Alpha = my_atof(ALPHA)/pow(10.0,my_atof(ALPHAp));
	  f.Beta = my_atof(BETA)/pow(10.0,my_atof(BETAp));
	  f.Gamma = my_atof(GAMMA)/pow(10.0,my_atof(GAMMAp));
	  ncon = (int) my_atoi(NCON);

	  remark_chars = (atoi(NRFAC)+atoi(NREM)+atoi(NDIS)+atoi(NERR));
	  if (remark_chars > 0)
	    remark_lines = remark_chars/80 + 1;
	  if (remark_chars % 80 == 0)
	    remark_lines = remark_chars/80;
	  
	  symmetry_lines = atoi(NOPR)/5 + 1;
#ifdef PRINT_HEADER
	  printf("A = %f  B = %f  C = %f  \n",f.A,f.B,f.C);
	  printf("NCON = %s\n",NCON);
	  printf("Alpha = %f Beta = %f Gamma = %f \n",f.Alpha,f.Beta,f.Gamma);
	  printf("remark_chars = %d remark_lines = %d symmetry_lines = %d \n",
		 remark_chars,remark_lines,symmetry_lines);
#endif
	  fill_orth_matrix(&f,&m);
	  
	  for (i = 1; i < 2; i++)
	    fgets(csd_line,sizeof(csd_line),file1);
	  for (i = 1; i <= symmetry_lines; i++)
	  {
	    fgets(csd_line,sizeof(csd_line),file1);
	  }
	  for (i = 1; i <= remark_lines; i++)
	  {
	    fgets(csd_line,sizeof(csd_line),file1);
	  }
	  k = 1;
	  for (i = 1; i <= Atoms; i += 3)
	  {
	    UpdateProgress();
	    fgets(csd_line,sizeof(csd_line),file1);
	    tokens = strlen(csd_line)/26;
	    for (j = 0; (j < tokens) && (k <= Atoms); j++)
	    {
	      start = j * 27;
	      my_strncpy(Type(k),&csd_line[start],5);
	      my_strncpy(Xstr,&csd_line[start + 5],7);
	      my_strncpy(Ystr,&csd_line[start + 12],7);
	      my_strncpy(Zstr,&csd_line[start + 19],7);
	      X(k) = my_atof(Xstr)/100000.0;
	      Y(k) = my_atof(Ystr)/100000.0;
	      Z(k) = my_atof(Zstr)/100000.0;
	      clean_atom_type(Type(k));
#ifdef DEBUG
	      printf("%d %s %f %f %f\n",k,Type(k),X(k),Y(k),Z(k));     
#endif
	      fract_to_cart(&Point(k),&m);
	    k++;
	    }
	  }
	}
      }

      if ((Atoms > 0) && (pending))
      {
	assign_radii(mol);
	assign_bonds(mol);
	assign_types(mol);
	build_connection_table(mol);
	assign_bond_order(mol);
      }
      if (!pending)
	strcpy(Title,REFCODE);
      pending = FALSE;
    }
  }

  /* find the start of the next record */
  pos = ftell(file1);
  while (fgets(csd_line,sizeof(csd_line), file1))
  {
    if (csd_line[0] == '#')
      break;
    pos = ftell(file1);
  }
  fseek(file1,pos,0);

  return(TRUE);
}


void my_strncpy(char *str1, char *str2, int len)
{
  strncpy(str1,str2,len);
  str1[len] = '\0';
}

double my_atof(char *the_str)
{
  if (strlen(the_str) > 0)
    return(atof(the_str));
  else 
    return(0.0);
}

double my_atoi(char *the_str)
{
  if (strlen(the_str) > 0)
    return(atoi(the_str));
  else 
    return(0);
}

int check_refcode(char *keywords, char *csd_line,ums_type *mol)
{
  int want_this = FALSE;
  char REFCODE[9];
  char the_keywords[BUFF_SIZE];
  
  strcpy(the_keywords,keywords);
  uppercase(the_keywords);
  my_strncpy(REFCODE,&csd_line[1],8);

  strcpy(REFCODE,gettoken(REFCODE," ",1));
  if ((!UseMenus) &&
      (EQ(the_keywords,"KEYWORDS GO HERE") ||
       (strstr(the_keywords,REFCODE) != NULL)))
    want_this = TRUE;
  if ((UseMenus) &&
      (EQ(the_keywords,"KEYWORDS GO HERE") ||
       (EQ(the_keywords,REFCODE))))
    want_this = TRUE;
  if ((strstr(the_keywords,"ALL")) || (strstr(the_keywords,"LAST")))
    want_this = TRUE;
/*  printf("%s %s %d\n",the_keywords,REFCODE,want_this); */
  return(want_this);
}



