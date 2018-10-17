/* $Id: rdg96.c,v 1.1 1996/06/25 12:25:29 wscott Exp $

   contributed by Walter Scott (wscott@igc.phys.chem.ethz.ch)

   This is a small routine to read a GROMOS96 formatted
   "position coordinate block" (POSITION) or a
   "reduced position coordinate block" (POSITIONRED)
   The former has name information (atom and residue names) while
   the latter only has coordinates.
   This version does not support the reading of binary
   GROMOS files.


   NOTE: Many programs specify the units of the coordinates (e.g. Angstrom).
   GROMOS96 does NOT, as all physical constants, from K_B to EPS are 
   NOT hardwired into the code, but specified by the user.
   This allows some (mostly Americans) to use GROMOS96 in KCal and
   Angstrom and the rest of us to use kJoule and nm.
   It also makes it easy to use reduced units.

   We get around this by supplying a routine, rd_sco_gr96, which
   will scale the coordinates by a factor after reading.
   This routine is then called with the factor set to 1.0 in 
   read_gr96A, or to 10.0 in read_gr96N depending on the users choice.
   Thus, the coordinates are always returned to Babel in Angstroms.
*/

#include "bbltyp.h"


/* read next noncomment line on file f in the_line
   return nonzero if successfull
*/
int gline(char* the_line, int nsz, FILE* f){
   char* res;
   res = fgets(the_line,nsz,f);
   while (res != 0 && the_line[0] == '#'){
      res = fgets(the_line,nsz,f);
   }
   return (res !=0);
}

/* return the number of noncomment lines before getting
   the END marker. rewind the file to the original position
*/
int cnt_atms(FILE* f){
   int atnum = 0;
   int pos,res,no_end;
   char the_line[BUFF_SIZE];
   char* endstr = "END\n";
/*begin*/
   pos =  ftell(f);
   res = gline(the_line,BUFF_SIZE,f);
   no_end = (strcmp(the_line,endstr) != 0);
   while (res && no_end){
      atnum++;
      res = gline(the_line,BUFF_SIZE,f);
      no_end = (strcmp(the_line,endstr) != 0);
   }
   fseek(f,pos,0);
   return atnum;
}

int rd_sco_gr96(FILE *file1, ums_type *mol,double fac){
  char the_line[BUFF_SIZE];
  int i,res,no_end,dummy,lfirst;
  char* endstr = "END\n";
  char* titstr = "TITLE\n";
  char* posstr = "POSITION\n";
  char* posredstr = "POSITIONRED\n";
  char* prgstr = "rd_sco_gr96";
  double x,y,z;

/*begin*/
  res = gline(the_line,BUFF_SIZE,file1);
  if (res ==0){
     fprintf(stderr,"%s: unexpected end of file\n",prgstr);
     return (FALSE);
  }

  Title[0] = 0;

  if (strcmp(the_line,titstr) ==0){
     fprintf(stderr,"got TITLE!!\n");fflush(stderr);

     /* we have a title block*/
     res = gline(the_line,BUFF_SIZE,file1);
     no_end = (strcmp(the_line,endstr) != 0);
     lfirst = TRUE;
     while (res && no_end){
	/* copy only the first line of the title block*/
	if (lfirst){
	   lfirst = FALSE;
	   strncpy(Title,the_line,200);
	}
	res = gline(the_line,BUFF_SIZE,file1);
	no_end = (strcmp(the_line,endstr) != 0);
     }
     res = gline(the_line,BUFF_SIZE,file1);
  }

/* do I really need to do this here ? */
  Bonds = 0;

  if (strcmp(the_line,posstr) ==0){
     /* we have a position block*/
     fprintf(stderr,"POSITION!!\n");fflush(stderr);
     Atoms = cnt_atms(file1);
     fprintf(stderr,"BLABLA %d\n",Atoms);fflush(stderr);
     initialize_ums(&mol);
     initialize_residues(&mol);
     ShowProgress(Atoms,"Reading Atoms");

     res = gline(the_line,BUFF_SIZE,file1);
     no_end = (strcmp(the_line,endstr) != 0);
     i = 1;
     while (res && no_end){
	res = sscanf(the_line,"%d %s %s %d %lf %lf %lf\n",
		     &ResNum(i),ResName(i),AtmId(i),&dummy,&x,&y,&z);
	if (res !=7){
	   fprintf(stderr,"%s: format error 2 (%d)\n",prgstr,res);
	   fprintf(stderr,"got '%s'\n",the_line);
	   return (FALSE);
	}
	X(i) = x*fac;
	Y(i) = y*fac;
	Z(i) = z*fac;

	i++;
	res = gline(the_line,BUFF_SIZE,file1);
	no_end = (strcmp(the_line,endstr) != 0);
     }
     if (no_end){
	fprintf(stderr,"%s: END expected\n",prgstr);
	return (FALSE);
     }else{
	fprintf(stderr,"OK %d\n",i);fflush(stderr);
	return (TRUE);
     }
  }else{
     if (strcmp(the_line,posredstr) ==0){
	fprintf(stderr,"POSRED!!\n");fflush(stderr);
	/* we have a positionred block*/
	Atoms = cnt_atms(file1);
	fprintf(stderr,"getting %d\n",Atoms);fflush(stderr);
	initialize_ums(&mol);
	initialize_residues(&mol);
	ShowProgress(Atoms,"Reading Atoms");


	res = gline(the_line,BUFF_SIZE,file1);
	no_end = (strcmp(the_line,endstr) != 0);
	i = 1;
	while (res && no_end){
	   res = sscanf(the_line,"%lf %lf %lf",&x,&y,&z);
	   if (res !=3){
	      fprintf(stderr,"%s: format error 1\n",prgstr);
	      return (FALSE);
	   }
	   X(i) = x*fac;
	   Y(i) = y*fac;
	   Z(i) = z*fac;

	   i++;
	   res = gline(the_line,BUFF_SIZE,file1);
	   no_end = (strcmp(the_line,endstr) != 0);
	}
	if (no_end){
	   fprintf(stderr,"%s: END expected\n",prgstr);
	   return (FALSE);
	}else{
	   return (TRUE);
	}
     }
  }
/*  assign_radii(mol);
  assign_bonds(mol); 
  assign_types(mol);
  build_connection_table(mol); 
  assign_bond_order(mol);
*/
  return(TRUE);
}

int read_gr96A(FILE *file1, ums_type *mol){
   return  rd_sco_gr96(file1,mol,1.0);
}


int read_gr96N(FILE *file1, ums_type *mol){
   return  rd_sco_gr96(file1,mol,10.0);
}

