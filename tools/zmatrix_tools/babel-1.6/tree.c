#include "bbltyp.h"

#define NewType(x)             new_mol->atoms[x].type
#define NewValence(x)          new_mol->atoms[x].valence
#define NewConnection(x,y)     new_mol->atoms[x].connected_atoms[y]
#define NewAtoms               new_mol->num_atoms
#define NewStart(x)            new_mol->connections[x].start
#define NewEnd(x)              new_mol->connections[x].end
#define NewBond_order(x)       new_mol->connections[x].bond_order
#define NewBonds               new_mol->num_bonds
#define NewX(a)                new_mol->atoms[a].point.x
#define NewY(a)                new_mol->atoms[a].point.y
#define NewZ(a)                new_mol->atoms[a].point.z

#define NewChainNum(x) new_mol->residues[x].chain_num
#define NewResNum(x) new_mol->residues[x].res_num
#define NewResName(x) new_mol->residues[x].res_type
#define NewAtmId(x) new_mol->residues[x].atm_type

static int cycle;

ums_type *renum_for_zmat(ums_type *mol,int base)
{
  z_tree *tree;
  cycle = 0;
  
  tree = (z_tree *)malloc((Atoms + 1) * sizeof(z_tree));
  build_z_tree(mol,base,tree);
  find_z_kids(mol,tree);
  dfs_z_tree(base,mol,tree);
  continuity_check(mol);
  mol = build_new_ums(mol,Atoms);
  push_hydrogens_to_end(mol);
  mol = build_new_ums(mol,Atoms);
  return(mol);
}

void build_z_tree(ums_type *mol, int root, z_tree *tree)
{
  int i,j,k;
  

  for (i = 1; i <= Atoms; i++)
  {
    Redo(i) = 0;
    tree[i].level = 0;
    tree[i].ancestor = 0;
  }
  
  Redo(root) = 1;
  
  for (i = 1; i <= Atoms; i++)
  {
    for (j = 1; j <= Atoms; j++)
      if (Redo(j) == i)
      {
	for (k = 0; k < Valence(j); k++)
	  if (Redo(Connection(j,k)) == 0)
	  {
	    (Redo(Connection(j,k)) = i + 1);
	    tree[Connection(j,k)].level = i + 1;
	    tree[Connection(j,k)].ancestor = j;
	  }
      }
  }
}  

void push_hydrogens_to_end(ums_type *mol)
{
  int i, num = 1;
  
  for (i = 1; i <= Atoms; i++)
    Redo(i) = 0;

  for (i = 1; i <= Atoms; i++)
  {
    if (Atomic_number(i) != 1)
    {
      Redo(i) = num;
      num++;
    }
  }

  for (i = 1; i <= Atoms; i++)
  {
    if (Redo(i) == 0)
    {
      Redo(i) = num;
      num++;
    }
  }
}

    

void continuity_check(ums_type *mol)
{
  int i;
  
  for (i = 1; i <= Atoms; i++)
  {
    if (Redo(i) == 0)
    {
      printf("%s %d\n",Type(i),i);
      cycle++;
      Redo(i) = cycle;
    }
  }
}

void find_z_kids(ums_type *mol,z_tree *tree)
{
  int i,j,k;
  
  for (i = 1; i <= Atoms; i++)
  {
    tree[i].kids = 0;
  }
  for (i = 1; i <= Atoms; i++)
  {
    for (j = 1; j <= Atoms; j++)
      if (tree[j].ancestor == i)
      {
	k = tree[i].kids;
	tree[i].kid[k] = j;
	tree[i].kids++;
      }
  }
}

void print_z_tree(ums_type *mol,z_tree *tree)
{
  int do_return;
  int i,j;
  
  for (i = 0; i <= Atoms; i++)
  {
    do_return = FALSE;
    for (j = 1; j <= Atoms; j++)
    {
      if (tree[j].level == i)
      {
	printf("%4d (%d) ",j,tree[j].ancestor);
	do_return = TRUE;
      }
    }
    if (do_return == TRUE)
      printf("\n");
  }

  for (i = 1; i <= Atoms; i++)
  {
    printf("%4d %4s- ",i,Type(i));
    for (j = 0; j < tree[i].kids; j++)
      printf("%4d",tree[i].kid[j]);
    printf("\n");
  }
}

void dfs_z_tree(int x, ums_type *mol, z_tree *tree)
{
  int i;
  cycle++;
  {
    Redo(x) = cycle;
    for (i = 0; i < tree[x].kids; i++)
    {
      if (tree[x].kids > 0)
	dfs_z_tree(tree[x].kid[i],mol,tree);
    }
  }
}

ums_type *renumber_ums(ums_type *mol,int heavy_count)
{
  int i;
  int j = 0;
  ums_type *new_mol;
  
  new_mol = (ums_type *)malloc(sizeof(ums_type));
  
  NewBonds = 0;
  NewAtoms = heavy_count;
  
  for (i = 0;i < Bonds; i ++)
  {
    if (Redo(Start(i)) != 0 && Redo(End(i)) != 0)
    {
      NewBonds ++;
    }
  }

  initialize_ums(&new_mol);
  if (HasResidues)
    initialize_residues(&new_mol);

  for (i = 1;i <= Atoms;i ++)
  {
    if (Redo(i) != 0)
    {     
      NewValence(Redo(i)) = 0;
      NewX(Redo(i)) = X(i);
      NewY(Redo(i)) = Y(i);
      NewZ(Redo(i)) = Z(i);
      strcpy(NewType(Redo(i)),Type(i));
      if (HasResidues)
      {
	NewChainNum(Redo(i)) = ChainNum(i);
	NewResNum(Redo(i)) = ResNum(i);
	strcpy(NewResName(Redo(i)),ResName(i));
	strcpy(NewAtmId(Redo(i)),AtmId(i));
      }
    }      
  }

  for (i = 0;i < Bonds; i ++)
  {
    if (Redo(Start(i)) != 0 && Redo(End(i)) != 0)
    {
      NewStart(j) = Redo(Start(i));
      NewEnd(j)   = Redo(End(i));
      j++;
    }
  }

  new_mol->control = mol->control;
  if (HasResidues)
    free(mol->residues);
  dissect_connection_table(new_mol);
  return(new_mol);
}









