#include "bbltyp.h"

#define DELIMS "\t\n "

int read_smiles(FILE *fp, ums_type *mol)
{
  smilescontab_t *ct = NULL;
  char smiles_string[BUFF_SIZE];
  char the_line[BUFF_SIZE];

  fgets(the_line,sizeof(the_line),fp);
  get_token(smiles_string,the_line,DELIMS,1);

  if (smiles_string)
    ct = smilestocontab(smiles_string);
  
  if (ct)
  {
    ct_to_ums(ct,mol);
    strcpy(Title,smiles_string);
    free(ct);
  }
  else
    return(FALSE);
  
  return(TRUE);
}

void ct_to_ums(smilescontab_t *ct, ums_type *mol)
{
  int i,j, atm;

  Atoms = ct->natoms;
  initialize_ums(&mol);
  
  for (i = 0; i < ct->natoms; i++)
  {
    atm = i+1;
    strcpy(Type(atm),ct->atom[i].symbol);
    for (j = 0; j < md_MAXBONDS; j++)
    {
      if (ct->atom[i].bondtype[j] == 0)
	break;
      else
      {
	Valence(atm)++;
	Connection(atm,j) = ct->atom[i].bondedto[j] + 1;
	BO(atm,j) = ct->atom[i].bondtype[j];
      }
    }
  }
  build_connection_table(mol);
  assign_type_by_bo(mol);
}

#undef DELIMS
