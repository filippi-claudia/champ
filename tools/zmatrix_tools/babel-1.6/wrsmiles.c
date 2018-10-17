#include "bbltyp.h"

int write_smiles(FILE *fp, ums_type *mol)
{
  smilescontab_t *ct;
  char *smilesstring;
  int result;
  
  ct = (smilescontab_t *)malloc(sizeof(smilescontab_t));
  
  find_aromatic_atoms(mol);
  mol = delete_atoms(mol,"default");
  if (ums_to_smiles_ct(mol,&ct))
  {
    smilesstring = contabtosmiles(ct);
    fprintf(fp,"%s\n",smilesstring);
  }
  else
  {
    printf("Unable to create smiles string\n");
    result = FALSE;
  }
  
  if (ct)
    free(ct);
  if (smilesstring)
    free(smilesstring);

  return(result);
}

/*-----------------------------------------------------
FUNCTION : ums_to_smiles_ct 
PURPOSE  : translate a ums into a SMILES connection table
------------------------------------------------------*/

int ums_to_smiles_ct(ums_type *mol, smilescontab_t **ct)
{
  smilescontab_t *newptr;
  int i,j;
  int conn, bo;
  
  if(!(newptr = realloc((*ct), smiles_CONTABHDRSIZE + (Atoms * sizeof(smilesatom_t)))))
  {
    printf("OUT OF MEMORY AT %s %s \n",__FILE__,__LINE__);
    return FALSE;
  }
  
  (*ct) = newptr;
  (*ct)->natoms = Atoms;

  for (i = 1; i <= Atoms; i++)
  {
    atomic_number_to_name(Atomic_number(i),(*ct)->atom[i-1].symbol);

    for (j = 0; j < md_MAXBONDS; j++)
    {
      (*ct)->atom[i-1].bondedto[j] = (j  < Valence(i) ? Connection(i,j) - 1 :  0);
      (*ct)->atom[i-1].bondtype[j] = (j  < Valence(i) ? BO(i,j) : 0);
    }
  }
  return(TRUE);
}


void printcontab(smilescontab_t *contab)
{
   int a, b;

   for(a=0; a<contab->natoms; a++)
   {
      printf("%4d %2.2s", a+1, contab->atom[a].symbol);
      for(b=0; b<md_MAXBONDS; b++)
         printf(" %4d (%d)", contab->atom[a].bondedto[b]+1, contab->atom[a].bondtype[b]);
      printf("\n");
   }
}
 






