/*
   smilesto.c

   converts SMILES strings to connection tables

   Simon Kilvington, University of Southampton, 1995
*/

#include "bbltyp.h"

#define isatomchr(C)		(isalpha(C) || (C) == '*')

/*
   smilestocontab
   uses a finite state machine to convert a SMILES string to a connection table
   returns NULL if any trouble (invalid SMILES or no memory) or a ptr to the alloc'd fragment
   only the first molecule specified in the string will be constructed, terminators are fullstop, space, newline, and end of string
   differences between the SMILES used by this and "real" SMILES strings are:
   1. charges, isotope numbers, and attached hydrogens can not be specified inside an atoms brackets, eg [NH4+]
   2. as coords are not generated symmetry symbols (ie @, @@, /, and \) can not be used
!!!!!!!!!!!! change it so these are just ignored rather than causing parsing errors !!!!!!!!!!!!!!!!!
*/

#define MAXSMILESRINGS 100	/* max ring closure number +1 */

smilescontab_t *
smilestocontab(char *smiles)
{
   int i, b, rno, loss, lastatom, maxbranches, bstackptr, *bstack = NULL, ringatom[MAXSMILESRINGS];
   char atname[2], *bord = NULL;
   smilescontab_t *frag = NULL;
   smilesatom_t *aptr = NULL;
   smilesbond_t nextbond;
   int toggle;

   static const char bondchr[] = "~-=#:";

/* the states */
   enum state_typ { St_Start, St_Name1, St_OpenSquare, St_Name2, St_CloseSquare, St_Number, St_CloseBranch, St_OpenBranch, St_Bond, St_End, St_Error } state;

/* the tokens */
   enum {Tok_Symbol, Tok_Bond, Tok_Number, Tok_OpenSquare, Tok_CloseSquare, Tok_OpenBranch, Tok_CloseBranch, Tok_End };

/* the state machine */
   static const int next[9][8] = { {St_Name1, St_Error, St_Error,  St_OpenSquare, St_Error,       St_Error,      St_Error,       St_End},
				   {St_Name1, St_Bond,  St_Number, St_OpenSquare, St_Error,       St_OpenBranch, St_CloseBranch, St_End},
				   {St_Name2, St_Error, St_Error,  St_Error,      St_Error,       St_Error,      St_Error,       St_Error},
				   {St_Error, St_Error, St_Error,  St_Error,      St_CloseSquare, St_Error,      St_Error,       St_Error},
				   {St_Name1, St_Bond,  St_Number, St_OpenSquare, St_Error,       St_OpenBranch, St_CloseBranch, St_End},
				   {St_Name1, St_Bond,  St_Number, St_OpenSquare, St_Error,       St_OpenBranch, St_CloseBranch, St_End},
				   {St_Name1, St_Bond,  St_Error,  St_OpenSquare, St_Error,       St_OpenBranch, St_CloseBranch, St_End},
				   {St_Name1, St_Bond,  St_Error,  St_OpenSquare, St_Error,       St_Error,      St_Error,       St_Error},
				   {St_Name1, St_Error, St_Error,  St_OpenSquare, St_Error,       St_Error,      St_Error,       St_Error}  };

   for(i=0; i<MAXSMILESRINGS; i++)
      ringatom[i] = -1;

   if(!(frag = malloc(smiles_CONTABHDRSIZE)))
   {
      printf("No memory to convert SMILES string\n");
      return NULL;
   }
   frag->natoms = 0;

   maxbranches = strutils_noccurrences(smiles, '(');
   if (maxbranches > 0)
     if(!(bstack = malloc(maxbranches * sizeof(int))))
     {
       printf("No memory to convert SMILES string\n");
       free(frag);
       return NULL;
     }
   
   i = 0;			/* index into smiles string */
   bstackptr = 0;		/* no of branches on the stack */
   state = St_Start;
   while(state != St_Error && state != St_End)
   {
   /* see what the next token in the string is, and move onto the appropriate state */
      if(strchr(". \n", smiles[i]) != NULL)		state = (enum state_typ) next[state][Tok_End];
      else if(isatomchr(smiles[i]))			state = (enum state_typ) next[state][Tok_Symbol];
      else if(bord = strchr(bondchr, smiles[i]))	state = (enum state_typ) next[state][Tok_Bond];
      else if(isdigit(smiles[i]) || smiles[i] == '%')	state = (enum state_typ) next[state][Tok_Number];
      else if(smiles[i] == '[')				state = (enum state_typ) next[state][Tok_OpenSquare];
      else if(smiles[i] == ']')				state = (enum state_typ) next[state][Tok_CloseSquare];
      else if(smiles[i] == '(')				state = (enum state_typ) next[state][Tok_OpenBranch];
      else if(smiles[i] == ')')				state = (enum state_typ) next[state][Tok_CloseBranch];
      else 						state = St_Error;			/* invalid character */

   /* do we need to do anything in this state */
      switch(state)
      {
      case St_Name1:
      /* check it's not "Cl" or "Br" */
	 if(strncmp(&smiles[i], "Cl", 2) == 0 || strncmp(&smiles[i], "Br", 2) == 0)
	 {
	    atname[0] = smiles[i];
	    atname[1] = smiles[i+1];
	    i+=2;
	 }
	 else
	 {
	    atname[0] = smiles[i];
            atname[1] = '\0';
	    i++;
	 }
	 if(!addatomtocontab(&frag, atname, &lastatom, nextbond))
	    state = St_Error;
	 nextbond = SMILESBOND_Single;
	 break;

      case St_Name2:
	 if(smiles[i+1] == ']')
	 {
	    atname[0] = smiles[i];
	    atname[1] = '\0';
	    i++;
	 }
	 else
	 {
	    atname[0] = smiles[i];
	    atname[1] = smiles[i+1];
	    i += 2;
	 }
	 if(!addatomtocontab(&frag, atname, &lastatom, nextbond))
	    state = St_Error;
	 nextbond = SMILESBOND_Single;
	 break;

      case St_Number:
	 if(smiles[i] == '%')
	 {
	    rno = (10 * (smiles[i+1] - '0')) + (smiles[i+2] - '0');
	    i += 3;
	 }
	 else
	 {
	    rno = smiles[i] - '0';
	    i++;
	 }
	 if(ringatom[rno] == -1)		/* start of ring */
	 {
	    ringatom[rno] = lastatom;
	 }
	 else					/* end of ring */
	 {
	    if(!closecontabring(frag, ringatom[rno], lastatom))
	       state = St_Error;
	    ringatom[rno] = -1;
	 }
	 break;

      case St_Bond:
	 nextbond = (bord - bondchr);
	 i++;
	 break;

      case St_OpenBranch:
	 bstack[bstackptr++] = lastatom;
	 i++;
	 break;

      case St_CloseBranch:
	 if(bstackptr == 0)
	 {
	    printf("Unmatched closing branch bracket in SMILES string at character %d\n", i+1);
	    state = St_Error;
	 }
	 else
	 {
	    lastatom = bstack[--bstackptr];
	    i++;
	 }
	 break;

      case St_End:
	 if(frag->natoms == 0)
	 {
	    printf("SMILES string contains no atoms\n");
	    state = St_Error;
	 }
	 break;

      case St_Error:
         printf("Error parsing SMILES string at character %d \"%c\"\n", i+1, smiles[i]);
	 break;

      default:
	 i++;
	 break;
      }
   }

/* check we are not in the middle of any rings or branches */
   if(state != St_Error && bstackptr > 0)
   {
      printf("Unmatched brackets in SMILES string\n");
      state = St_Error;
   }

   for(i=0; state != St_Error && i<MAXSMILESRINGS; i++)
   {
      if(ringatom[i] != -1)
      {
	 printf("Unmatched ring closure number (%d) in SMILES string\n", i);
	 state = St_Error;
      }
   }

   if(state == St_Error)
   {
   /* clean up */
     if (frag)
     {
       free(frag);
       frag = NULL;
     }
     
   }
   else
   {
   /* convert all the chemical symbols to the correct case */
      for(i=0; i<frag->natoms; i++)
      {
	 aptr = &frag->atom[i];
	 aptr->symbol[0] = toupper(aptr->symbol[0]);
	 aptr->symbol[1] = tolower(aptr->symbol[1]);
      }
   }

   if (bstack)
   {
     free(bstack);
     bstack = NULL;
   }

   return frag;
}

#undef MAXSMILESRINGS

/*
   addatomtocontab
   adds an extra atom to *frag and bonds it to *lastatom (unless the new atom is the first one) with the given bond type
   if bondtype is SINGLE and both atoms have lower case names an AROMATIC bond is used instead
   updates *lastatom to be the index of the atom just added (ie the last atom in the fragment)
   returns FALSE if any trouble (no memory, unknown atom type, or too many bonds to *lastatom)
*/

int
addatomtocontab(smilescontab_t **frag, char *atname, int *lastatom, smilesbond_t bondtype)
{
   smilescontab_t *newptr;
   int type, newatom, b, btype, bbits;
   smilesatom_t *aptr, *bptr;
   char upper[2], pdbatm[8];
   int uppercase, err;

   err = FALSE;

   uppercase = !aromaticsmilessym(atname);
   upper[0] = toupper(atname[0]);
   upper[1] = toupper(atname[1]);

   if(!uppercase && bondtype == SMILESBOND_Single)
   {
      if(aromaticsmilessym((*frag)->atom[*lastatom].symbol))
         bondtype = SMILESBOND_Aromatic;
   }

   newatom = (*frag)->natoms;

   if(!(newptr = realloc((*frag), smiles_CONTABHDRSIZE + ((*frag)->natoms + 1) * sizeof(smilesatom_t))))
   {
      printf("No memory to build fragment\n");
      return FALSE;
   }
   (*frag) = newptr;
   (*frag)->natoms ++;

/* set up "newatom" and bond it to "*lastatom" (unless newatom is the first) */
   aptr = &((*frag)->atom[newatom]);
   strncpy(aptr->symbol, atname, 2);
   if(newatom == 0)
   {
      for(b=0; b<md_MAXBONDS; b++)
      {
	 aptr->bondedto[b] = md_NOBOND;
	 aptr->bondtype[b] = SMILESBOND_NoBond;
      }
   }
   else
   {
      bptr = &((*frag)->atom[*lastatom]);
      b = nextfreebondto(bptr);
      if(b != -1)
      {
	 bptr->bondedto[b] = newatom;
	 bptr->bondtype[b] = bondtype;
	 aptr->bondedto[0] = *lastatom;
	 aptr->bondtype[0] = bondtype;
	 for(b=1; b<md_MAXBONDS; b++)
	 {
	    aptr->bondedto[b] = md_NOBOND;
	    aptr->bondtype[b] = SMILESBOND_NoBond;
	 }
      }
      else
      {
	 printf("Too many bonds to atom %d\n", (*lastatom)+1);
	 err = TRUE;
      }
   }

/* set up the return values */
   *lastatom = newatom;

   return !err;
}

/*
   closecontabring
   bonds atom1 and atom2 using either a SINGLE or AROMATIC bond depending on whether the names of both are lower case or not
   returns FALSE if any trouble (too many bonds to one of the atoms)
*/

int
closecontabring(smilescontab_t *frag, int atom1, int atom2)
{
   int bto1, bto2;
   smilesatom_t *a1ptr, *a2ptr;
   smilesbond_t btype;
   int lower1, lower2;

   a1ptr = &frag->atom[atom1];
   a2ptr = &frag->atom[atom2];

   if((bto1 = nextfreebondto(a1ptr)) == -1)
   {
      printf("Too many bonds to atom %d\n", atom1+1);
      return FALSE;
   }
   if((bto2 = nextfreebondto(a2ptr)) == -1)
   {
      printf("Too many bonds to atom %d\n", atom2+1);
      return FALSE;
   }

   lower1 = aromaticsmilessym(a1ptr->symbol);
   lower2 = aromaticsmilessym(a2ptr->symbol);

   btype = (lower1 && lower2) ? SMILESBOND_Aromatic : SMILESBOND_Single;

   a1ptr->bondedto[bto1] = atom2;
   a1ptr->bondtype[bto1] = btype;

   a2ptr->bondedto[bto2] = atom1;
   a2ptr->bondtype[bto2] = btype;

   return TRUE;
}

/*
   aromaticsmilessym
   returns TRUE if the given SMILES atom symbol is lower case, and therefore aromatic
*/

int
aromaticsmilessym(char *sym)
{
   return (sym[0] == tolower(sym[0]));		/* as metals may be called 'Zn' etc */
}

/*
   nextfreebondto
   returns index of first bondedto[] that is md_NOBOND, or -1 if there are no free bonds
*/

int
nextfreebondto(smilesatom_t *aptr)
{
   int f;
   int freebond = FALSE;

   for(f=0; f<md_MAXBONDS && !freebond; f++)
      freebond = (aptr->bondedto[f] == md_NOBOND);

   return (freebond) ? f-1 : -1;
}

/*
   strutils_noccurrances
   counts the number of times the given character appears in the given string
*/

int
strutils_noccurrences(char *buffer, char c)
{
   int i, n;

   i = n = 0;
   while(buffer[i] != '\0')
   {
      n += (buffer[i] == c);
      i++;
   }

   return n;
}
