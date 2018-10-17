/*
   tosmiles.c

   converts connection tables to unique SMILES strings

   Simon Kilvington, University of Southampton, 1995
*/

#include "bbltyp.h"

/*
   contabtosmiles
   converts the given fragment into a unique SMILES string
   uses the algorithm described by the Weiningers in J.Chem.Inf.Comput.Sci. 29 (1989), 97-101
   returns a ptr to the string which needs to be free'd when you've finished with it
   returns NULL if any problems
   it will explicitly include any H atoms in frag in the resulting SMILES, so it's best to remove all H atoms from frag before calling this
*/

char *
contabtosmiles(smilescontab_t *frag)
{
   char *smiles;
   int *rank, *primesum, *tmprank, *primelist;
   int a, b, c, n, bto, nranks, lowest;
   int nbonds, nprimes;
   int prime, same;
   smilesatom_t *aptr;
   block_ptr blk;
   int i;

   if(!block_alloc(&blk, "iiii", &rank, frag->natoms, &primesum, frag->natoms, &tmprank, frag->natoms, &primelist, frag->natoms))
      return NULL;

/* generate an array of the first "frag->natoms" prime numbers */
   nprimes = 0;
   for(n=2; nprimes<frag->natoms; n++)
   {
      prime = TRUE;
      for(c=2; c<=(n/2) && prime; c++)
	 prime = ((n % c) != 0);
      if(prime)
         primelist[nprimes++] = n;
   }

/* determine the initial graph invariants for each atom, store them in primesum */
   for(a=0; a<frag->natoms; a++)
   {
      aptr = &frag->atom[a];
   /* no of connections */
      nbonds = nbondsto(aptr);
      primesum[a] = nbonds * 100000;
   /* no of non-H bonds (frag should contain no H atoms) */
      primesum[a] += nbonds * 10000;
#if 0	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   /* atomic no */
      primesum[a] += md_atomicnumber(aptr->symbol) * 100;
   /* sign and absolute value of charge */
      primesum[a] += md_chargeonatom(aptr) * 10;
   /* no of attached H's */
      primesum[a] += (md_maxvalency(aptr->symbol) - nbonds);
#endif	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   }

   do
   {
   /* rank each atom based on primesum values (initially these are the invariants, on subsequent loops they are the tie breaker ranks) */
      sortranks(frag->natoms, primesum, rank);
      do
      {
      /* generate the product of the "corresponding primes" of each atoms neighbours */
	 for(a=0; a<frag->natoms; a++)
	 {
	    aptr = &frag->atom[a];
	    primesum[a] = 1;
	    for(b=0; b<md_MAXBONDS; b++)
	    {
	       bto = aptr->bondedto[b];
	       primesum[a] *= (bto != -1) ? primelist[rank[bto]-1] : 1;		/* -1 as array indices start at 0 */
	    }
	 }
   
      /* resolve ties in "rank" list using the "primesum" values, put the new ranks in "tmprank" */
	 nranks = resolveties(frag->natoms, rank, primesum, tmprank);
   
      /* see if the ranks are the same as last time */
	 same = TRUE;
	 for(a=0; a<frag->natoms && same; a++)
	    same = (rank[a] == tmprank[a]);

      /* move the sorted data back into rank array */
         memcpy(rank, tmprank, frag->natoms * sizeof(int));
      }
      while(!same);
   
   /* see if we have any tied atoms */
      if(nranks != frag->natoms)		/* break ties */
      {
      /* find lowest tied rank */
	 lowest = 0;
	 do
	 {
	    lowest++;
	    c = 0;
	    for(a=0; a<frag->natoms; a++)
	       c += (rank[a] == lowest);
	 }
	 while(c == 1);
      /* double all the ranks, store them in primesum array */
	 for(a=0; a<frag->natoms; a++)
	    primesum[a] = rank[a] * 2;
      /* reduce the first lowest tied one by 1 */
	 for(a=0; primesum[a] != (lowest * 2); a++)
	    ;
	 primesum[a]--;
      }

   }
   while(nranks != frag->natoms);

/* frag is now cannonicalised, lets get on with generating the smiles string */

   smiles = generatesmilesstring(frag, rank);

/* clean up */
   block_free(&blk);

   return smiles;
}

/*
   sortranks
   generates a list of ranks in "rank" based on the numbers in "value" (ie lowest value is ranked 1 etc)
   returns number of different ranks
   the value array will be corrupted, it'll end up as all -1's
*/

int
sortranks(int natoms, int *value, int *rank)
{
   int a, sort, lowest;

   sort = 1;
   do
   {
   /* find the lowest value that isnt -1, our tag for values that have already been sorted out (ho ho!) */
      lowest = -1;		/* not found yet */
      for(a=0; a<natoms; a++)
	 lowest = ((lowest == -1) || ((value[a] != -1) && (value[a] < lowest))) ? value[a] : lowest;
   /* pick out all the lowest values and sort them based on their previous ranks */
      if(lowest != -1)
      {
      /* all tied values have the same rank */
	 for(a=0; a<natoms; a++)
	 {
	    if(value[a] == lowest)
	    {
	       rank[a] = sort;
	       value[a] = -1;
	    }
	 }
	 sort++;
      }
   }
   while(lowest != -1);

   return sort-1;
}

/*
   resolveties
   resolves ties in "oldrank" array by ranking ties based on their "primesum" values
   new ranks are placed in "rank" array
   returns the new number of ranks
*/

int
resolveties(int natoms, int *oldrank, int *primesum, int *rank)
{
   int a, tie, lowest, nties, nranks;

/* initialise */
   nranks = 1;
   for(a=0; a<natoms; a++)
      rank[a] = -1;			/* this entry not ranked yet */

/* find the entry that has oldrank==tie and the lowest primesum, set its rank entry to nranks */
   for(tie=1; tie<=natoms; tie++)
   {
      do
      {
	 nties = 0;
	 lowest = -1;
	 for(a=0; a<natoms; a++)
	 {
	    if(rank[a] == -1 && oldrank[a] == tie)
	    {
	       nties++;
	       lowest = ((lowest==-1) || (primesum[a] < lowest)) ? primesum[a] : lowest;
	    }
	 }
      /* pick out all entries with oldrank==tie and primesum==lowest */
	 if(nties > 0)
	 {
	    for(a=0; a<natoms; a++)
	    {
	       if(oldrank[a] == tie && primesum[a] == lowest)
		  rank[a] = nranks;
	    }
	    nranks++;
	 }
      }
      while(nties > 0);
   }

   return nranks-1;
}

/*
   generatesmilesstring
   generates a smiles string for the fragment based on the given atomic ranks
*/

typedef struct				/* an array of these holds the data about each atom generated by the first pass thru the fragment */
{					/* the second pass uses these to write out the smiles string */
   char          symbol[2];		/* atom symbol */
   int          aromatic;		/* TRUE if ANY of the bonds to it are aromatic */
   int           ring[md_MAXBONDS];	/* ring closure numbers associated with this atom */
   int          rbracket;		/* does this atom end a branch */
   int          lbracket;		/* does the NEXT atom start a branch */
   int          endfrag;		/* do we need a "." after this atom to separate two structures */
   smilesbond_t bond;			/* order of bond to next atom */
} _smileyatomdata;

char *
generatesmilesstring(smilescontab_t *frag, int *rank)
{
   char *smiles;
   char tmpstr[8];
   int *bracketstack, *ringstack, *neworder;		/* neworder[i] is the the index of frag->atom[i] in the atomdata array */
   int bsptr, rsptr;
   int a, b, bto, f, lowest, tmp, smileslen;
   int smipos, atom, nextatom, ringatom, newring, closeatom, maxbonds, natoms, nrings;
   smilesbond_t order;
   int *visited;
   block_ptr blk;
   _smileyatomdata *atomdata, initatomdata;

   static const char bondchr[] = "-=#:";

   natoms = frag->natoms;

   if(!(atomdata = malloc(natoms * sizeof(_smileyatomdata))))
      return NULL;

   if(!block_alloc(&blk, "iiiB", &bracketstack, natoms, &ringstack, natoms*md_MAXBONDS, &neworder, natoms, &visited, natoms))
   {
      free(atomdata);
      return NULL;
   }

/* set up the initial (blank) atomdata entry */
   initatomdata.symbol[0] = '\0';
   initatomdata.symbol[1] = '\0';
   initatomdata.aromatic = FALSE;
   for(b=0; b<md_MAXBONDS; b++)
      initatomdata.ring[b] = -1;
   initatomdata.rbracket = FALSE;
   initatomdata.lbracket = FALSE;
   initatomdata.endfrag = FALSE;
   initatomdata.bond = SMILESBOND_NoBond;

   bsptr = 0;				/* points to next free entry on bracket stack */
   rsptr = 0;				/* ditto ring closure stack */
   for(a=0; a<frag->natoms; a++)	/* haven't visited any atoms yet */
      visited[a] = FALSE;
   natoms = 0;				/* number of atoms we have in our smiley atom data array so far */

/* find the atom with rank 1 */
   atom = lowestunvisitedatom(frag->natoms, rank, visited);

/* first pass, determine which order the atoms should appear in, and store some attributes of each one - start with the one we have just found */
   while(natoms < frag->natoms && atom != -1)
   {
   /* set up the atomdata entry for atom number "atom" */
      atomdata[natoms] = initatomdata;
      strncpy(atomdata[natoms].symbol, frag->atom[atom].symbol, 2);
      atomdata[natoms].aromatic = hasaromaticbonds(&frag->atom[atom]);
   /* remember we have visited it, and store its index into the atomdata array */
      visited[atom] = TRUE;
      neworder[atom] = natoms;
   /* if there are no unvisited atoms attached to this atom */
      if(nunvisitedbondsto(&frag->atom[atom], visited) == 0)
      {
      /* jump back to the atom on top of the branch stack if there is one */
	 if(bsptr > 0)
	 {
         /* add an end of branch bracket */
            atomdata[natoms].rbracket = TRUE;
         /* go back to the start of the branch */
	    bsptr --;
	    atom = bracketstack[bsptr];
	 }
         else
	 {
	 /* is it time to start a new structure, or have we reached the end of the molecule */
	    if(natoms != frag->natoms-1)
	    {
	       atomdata[natoms].endfrag = TRUE;
	       atom = lowestunvisitedatom(frag->natoms, rank, visited);
	    }
         /* increase the number of atoms we have seen, ready for next time */
            natoms ++;
	    continue;
	 }
      }
   /* there are unvisited atoms attached to this one */
/*!!!!!!!! what if the atom we pulled off the branch stack has no unvisited atoms attached to it?  !!!!!!!!!!*/
/*!!!!!!!! this can only happen if it's the last atom, or the end of a substructure                !!!!!!!!!!*/
nextatom=-1;
   /* find the lowest ranked, unvisited, atom attached to "atom", this is the next one to visit */
      lowest = frag->natoms+1;
      for(b=0; b<md_MAXBONDS; b++)
      {
	 bto = frag->atom[atom].bondedto[b];
	 if(bto != md_NOBOND && !visited[bto] && rank[bto] < lowest)
	 {
	    lowest = rank[bto];
	    nextatom = bto;
	 }
      }
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if(nextatom==-1)
{
   printf("natoms=%d\tatom=%d\tnext=-1\n", natoms, atom+1);
   atom=-1;
   continue;
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   /* count how many bonds from this atom form rings */
      nrings = 0;
   /* see if "nextatom" forms a ring, ie can we follow "nextatom" and get back to "atom" again */
      ringatom = findunvisitedring(frag, atom, nextatom, visited);
   /* did we find a ring */
      if(ringatom != -1)
      {
      /* is the bond order between "atom" and "ringatom" double or triple */
	 order = getbondorder(&frag->atom[atom], ringatom);
	 if(order == SMILESBOND_Double || order == SMILESBOND_Triple)
	 {
	 /* follow "ringatom" instead of "nextatom" so we dont have to close the ring with a multiple bond */
	    tmp = nextatom;
	    nextatom = ringatom;
	    ringatom = tmp;
	 }
      /* store the two atoms on the ring closure stack */
	 ringstack[rsptr++] = atom;
	 ringstack[rsptr++] = ringatom;
      /* increase our ring count */
	 nrings ++;
      /* see if any of the other unvisited atoms connected to this one also form a ring */
	 for(b=0; b<md_MAXBONDS; b++)
	 {
	 /* check there is an atom here, its not the ring we found before, and its unvisited */
	    bto = frag->atom[atom].bondedto[b];
	    if(bto == md_NOBOND || bto == nextatom || bto == ringatom || visited[bto])
	       continue;
	 /* does atom "bto" form a ring (ie can we trace a path from it back to "atom") */
	    newring = findunvisitedring(frag, atom, bto, visited);
	    if(newring != -1)
	    {
	    /* add the two atoms to the ring stack */
	       ringstack[rsptr++] = atom;
	       ringstack[rsptr++] = bto;
	    /* count how many bonds form rings */
	       nrings ++;
	    }
	 }
      }
   /* see if this is a branch point, ie more than 1 unvisited bonds to it that don't form rings */
      maxbonds = 1 + nrings;
      if(nunvisitedbondsto(&frag->atom[atom], visited) > maxbonds)
      {
	 atomdata[natoms].lbracket = TRUE;
      /* store this atom on the branch stack */
	 bracketstack[bsptr] = atom;
	 bsptr ++;
      }
   /* store the bond order between this atom and the next */
      atomdata[natoms].bond = getbondorder(&frag->atom[atom], nextatom);
   /* move onto the next atom */
      atom = nextatom;
   /* increase the number of atoms we have seen, ready for next time */
      natoms ++;
   }

/* did anything go wrong */
   if(natoms != frag->natoms)
   {
      free(atomdata);
      block_free(&blk);
      return NULL;
   }

/* pull the ring closure atoms off the ringstack and convert them into ring closure numbers */
   nrings = rsptr / 2;
   while(nrings > 0)
   {
      ringatom  = neworder[ringstack[--rsptr]];		/* convert the frag->atom index into an atomdata index */
      closeatom = neworder[ringstack[--rsptr]];
   /* mark the ring closure in atomdata[] */
      for(f=0; atomdata[ringatom].ring[f] != -1; f++)	/* find next free entry */
	    ;
      atomdata[ringatom].ring[f] = nrings;
      for(f=0; atomdata[closeatom].ring[f] != -1; f++)
	    ;
      atomdata[closeatom].ring[f] = nrings;
      nrings--;
   };

/* work out the max length "smiles" needs to be */
   smileslen = 1;					/* don't forget the terminator */
   for(a=0; a<frag->natoms; a++)
   {
   /* the atom name */
      smileslen += addsmilesatom(tmpstr, atomdata[a].symbol, atomdata[a].aromatic);
   /* ring closure numbers */
      for(b=md_MAXBONDS-1; b>=0; b--)
      {
         nrings = atomdata[a].ring[b];
	 if(nrings != -1)
	    smileslen += (nrings > 9) ? 3 : 1;		/* do we need a "%xx" or just a "x" */
      }
   /* right brackets */
      if(atomdata[a].rbracket)
	 smileslen ++;
   /* left brackets */
      if(atomdata[a].lbracket)
	 smileslen ++;
   /* bond symbols - we never need to write aromatic bond symbols because the bound atoms name's will be in lower case */
      order = atomdata[a].bond;
      if(order == SMILESBOND_Double || order == SMILESBOND_Triple)
	 smileslen ++;
   /* substructure separators */
      if(atomdata[a].endfrag)
         smileslen ++;
   }

/* alloc some space for "smiles" */
   if(!(smiles = malloc(smileslen)))
   {
      free(atomdata);
      block_free(&blk);
      return NULL;
   }

/* second pass - generate the smiles string using the data we have just set up */
   smipos = 0;				/* smipos is index of terminator in smiles string */
   for(a=0; a<frag->natoms; a++)
   {
   /* print the atom name */
      smipos += addsmilesatom(&smiles[smipos], atomdata[a].symbol, atomdata[a].aromatic);
   /* ring numbers are stored in reverse numerical order, so print them out backwards */
      for(b=md_MAXBONDS-1; b>=0; b--)
      {
         nrings = atomdata[a].ring[b];
	 if(nrings != -1)
	    smipos += addsmilesring(&smiles[smipos], nrings);
      }
   /* does this atom end a branch */
      if(atomdata[a].rbracket)
	 smiles[smipos++] = ')';
   /* does the atom end a substructure */
      if(atomdata[a].endfrag)
         smiles[smipos++] = '.';
   /* does the next atom start a branch */
      if(atomdata[a].lbracket)
	 smiles[smipos++] = '(';
   /* we never need to write aromatic bond symbols because the bound atoms name's will be in lower case */
      order = atomdata[a].bond;
      if(order == SMILESBOND_Double || order == SMILESBOND_Triple)
	 smiles[smipos++] = bondchr[order - SMILESBOND_Single];
   }
   smiles[smipos] = '\0';

   free(atomdata);
   block_free(&blk);

   return smiles;
}

/*
   addsmilesatom
   adds the smiles symbol(s) (in lower case if aromatic is TRUE) to the string
   returns the number of characters added
*/

#define NINSUBSET 10

int
addsmilesatom(char *string, char *symbol, int aromatic)
{
   int n, s;

   static const char *organicsubset[NINSUBSET] = { "B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I" };

/* if the symbol is not in the "organic subset" we need to enclose it in square brackets */
   for(s=0; s < NINSUBSET && strncmp(symbol, organicsubset[s], 2) != 0; s++)
      ;

/* count how many chars we add to the string */
   n = 0;

   if(s == NINSUBSET)
      string[n++] = '[';

   string[n++] = (aromatic) ? tolower(symbol[0]) : symbol[0];

/* is it a one letter symbol */
   if(symbol[1] != '\0')
      string[n++] = (aromatic) ? tolower(symbol[1]) : symbol[1];

   if(s == NINSUBSET)
      string[n++] = ']';

   return n;
}

#undef NINSUBSET

/*
   addsmilesring
   adds a ring closure digit (or %xx if ringnum > 9) to the given string
   returns the number of characters it added
*/

int
addsmilesring(char *string, int ringnum)
{
   int n;

   if(ringnum < 10)
   {
      *string = '0' + ringnum;
      n = 1;
   }
   else
   {
      n = (int)sprintf(string, "%%%d", ringnum);
   }

   return n;
}

/*
   lowestunvisitedatom
   returns the index of the atom which has "visited[]" set to FALSE and the lowest value in "rank[]"
*/

int
lowestunvisitedatom(int natoms, int *rank, int *visited)
{
   int a, lowestatom, lowestrank;

   lowestatom = -1;
   for(a=0; a<natoms; a++)
   {
      if(visited[a])
         continue;
      if(lowestatom == -1 || rank[a] < lowestrank)
      {
	 lowestrank = rank[a];
	 lowestatom = a;
      }
   }

   if(lowestatom == -1)
      printf("No unvisited atoms left\n");

   return lowestatom;
}

/*
   findunvisitedring
   determines whether atoms start and branch are in an unvisited ring, ie one with visited[i]==FALSE for all atoms i in the ring
   if they are it returns the other atom attached to start that is in the ring
   if they are not in a ring it returns -1
*/

int
findunvisitedring(smilescontab_t *frag, int start, int branch, int *visited)
{
   int b, close;
   int ring, *done;

/* get some space for a copy of the visited array as we will need to change it */
/*!!!!!! do we really need to copy it? what if "dotheymeet" restored "visited[bto]" after each call? !!!!!!*/
   if(!(done = malloc(frag->natoms * sizeof(int))))
      return -1;

   ring = FALSE;
   for(b=0; b<md_MAXBONDS && !ring; b++)
   {
      close = frag->atom[start].bondedto[b];
      if(close != md_NOBOND && close != branch && !visited[close])
      {
      /* take a copy of the visited array as we will need to change it */
         memcpy(done, visited, frag->natoms * sizeof(int));
      /* can we trace a path from "start", along the bond to "close", and end up at "branch" ? */
	 ring = dotheymeet(frag, start, close, branch, done);
      }
   }

   free(done);

   return ring ? close : -1;
}

/*
   dotheymeet
   recursively searches from "next" (but not along the bond to "start") for "target"
   only visits an atom i if visited[i]==FALSE
   all atoms that are visited will have visited[i] set to TRUE
   returns TRUE if it was found
*/

int
dotheymeet(smilescontab_t *frag, int start, int next, int target, int *visited)
{
   int b, bto;
   int meet;

/*!!!!! could we store "visited" for each atom connected to "next" here, then restore it before we return? !!!!!*/

   meet = FALSE;
   for(b=0; b<md_MAXBONDS && !meet; b++)
   {
      bto = frag->atom[next].bondedto[b];
      if(bto != md_NOBOND && bto != start && !visited[bto])
      {
	 visited[bto] = TRUE;
	 meet = (bto == target) ? TRUE : dotheymeet(frag, next, bto, target, visited);
      }
   }

   return meet;
}

/*
   nunvisitedbondsto
   returns the number of bonds to the given atom
   bonds are only counted if visited[i] == FALSE where i is the index of the atom bound to this one
*/

int
nunvisitedbondsto(smilesatom_t *aptr, int *visited)
{
   int b, n, bto;

   n=0;
   for(b=0; b<md_MAXBONDS; b++)
   {
      bto = aptr->bondedto[b];
      n += ((bto != md_NOBOND) && !visited[bto]);
   }

   return n;
}

/*
   nbondsto
   returns the total number of atoms bonded to this one
*/

int
nbondsto(smilesatom_t *aptr)
{
   int b, nbonds;

   nbonds = 0;
   for(b=0; b<md_MAXBONDS; b++)
      nbonds += (aptr->bondedto[b] != md_NOBOND);

   return nbonds;
}

/*
   getbondorder
   returns the bond order of the bond from the atom pointed to by "aptr" to atom number "atom"
   if they are not bonded it returns SMILESBOND_NoBond
*/

smilesbond_t
getbondorder(smilesatom_t *aptr, int atom)
{
   int b;
   smilesbond_t order;

   order = SMILESBOND_NoBond;
   for(b=0; b<md_MAXBONDS && order==SMILESBOND_NoBond; b++)
   {
      if(aptr->bondedto[b] == atom)
	 order = aptr->bondtype[b];
   }

   return order;
}

/*
   hasaromaticbonds
   returns TRUE if there are one or more aromatic bonds to the given atom
*/

int
hasaromaticbonds(smilesatom_t *aptr)
{
   int b;
   int has = FALSE;

   for(b=0; b<md_MAXBONDS && !has; b++)
      has = (aptr->bondedto[b] != md_NOBOND && aptr->bondtype[b] == SMILESBOND_Aromatic);

   return has;
}
