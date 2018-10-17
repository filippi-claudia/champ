/* ring detection routines written by Pat Walters */

#include "bbltyp.h"

int find_rings(ums_type *mol, ring_struct *rings)
{
  spanning_tree *tree1, *tree2, *tree3, *common_tree;
  connect_type *rca;
  int frj;
  int ring_count;
  int i,l;
  path path1,path2;
  ring_struct temp_list;
  set_type **set2, **set3, *common_set;
  int common;
  int path_size;

  frj =  Bonds - Atoms + 1;

  rings->ring_list = (path *)malloc(4 * Atoms * sizeof(path));
  rings->count = 0;
  rca = (connect_type *)malloc(frj * sizeof(connect_type));
  memset(rca,0,frj * sizeof(connect_type));
  tree1 = (spanning_tree *)malloc((Atoms+1) * sizeof(spanning_tree));
  tree2 = (spanning_tree *)malloc((Atoms+1) * sizeof(spanning_tree));
  tree3 = (spanning_tree *)malloc((Atoms+1) * sizeof(spanning_tree));
  common_tree = (spanning_tree *)malloc((Atoms) * sizeof(spanning_tree));
  common_set = init_set_minbits(Atoms);
  set2 = (set_type **)malloc((Atoms) * sizeof(set_type *));
  set3 = (set_type **)malloc((Atoms) * sizeof(set_type *));

  for (i = 0; i < Atoms; i++)
  {
    set2[i] = init_set_minbits(Atoms);
    set3[i] = init_set_minbits(Atoms);
  }

  path1.path_atoms = (int *)malloc((Atoms) * sizeof(int));
  path1.path_set = init_set_minbits(Atoms);
  path2.path_atoms = (int *)malloc((Atoms) * sizeof(int));
  path2.path_set = init_set_minbits(Atoms);

  temp_list.ring_list = (path *)malloc(Atoms * sizeof(path));
  memset(temp_list.ring_list,0,Atoms * sizeof(path));

  build_spanning_tree(mol,1,tree1);
  
  ring_count = find_closure_bonds(mol,tree1,rca);


#ifndef ULTRIX
  qsort(rca,ring_count,sizeof(connect_type),QSORT_PROTO sort_rca);
#else
  qsort(rca,ring_count,sizeof(connect_type), sort_rca);
#endif

  for (i = 0; i < ring_count; i++)
  {
    build_restricted_tree(mol,rca[i].start,rca[i].end,tree2);
    tree_to_sets(mol,tree2,set2);
    build_restricted_tree(mol,rca[i].end,rca[i].start,tree3);
    tree_to_sets(mol,tree3,set3);

    find_common_atom(mol,set2,set3,common_set);
    path_size = build_common_array(common_set,tree3,common_tree);
    common = 0;

    for (l = 0;l < path_size; l++)
    {
      common = common_tree[l].ancestor;
      init_path(Atoms,&path1);
      init_path(Atoms,&path2);
      path_to_root(tree2,common,&path1);
      path_to_root(tree3,common,&path2);
      paths_to_ring(path1,path2,&temp_list.ring_list[l],mol->num_atoms);
    }
    sort_rings(temp_list.ring_list,l);
    find_bogus_rings(temp_list.ring_list,l,mol->num_atoms); 
	save_good_rings(rings,&temp_list,l,TRUE);
  }

  sort_rings(rings->ring_list,rings->count);
  find_bogus_rings2(mol,rings->ring_list,rings->count,frj); 

  i =	rings->count;
  rings->count = 0;
  save_good_rings(rings,rings,i,TRUE);

  rings->ring_list = (path *)realloc(rings->ring_list,
				     (sizeof(path) * rings->count));
	
  if (rca) free(rca); 
  if (tree1) free(tree1); 
  if (tree2) free(tree2);
  if (tree3) free(tree3);
  if (common_tree) free(common_tree);
  free_set(common_set);

  for (i = 0;i < Atoms;i++)
  {
    free_set(set2[i]); 
    free_set(set3[i]);
  }
  free(set2);
  free(set3);

  if (path1.path_atoms) free(path1.path_atoms);
  if (path2.path_atoms) free(path2.path_atoms);
  free_set(path1.path_set);
  free_set(path2.path_set);
	
  free(temp_list.ring_list);

  if (rings->count != frj)
  {
    printf("ring count != frerejaque number\n");
    printf("Atoms = %d Bonds = %d ring count = %d frj = %d\n",
	   Atoms,Bonds,rings->count,frj);
    exit(0);
  }

  return(rings->count);
}

void tree_to_sets(ums_type *mol,spanning_tree *the_tree,set_type *the_set[])
{
  int i;
  
  for (i = 0;i < Atoms;i++)
    setclear(the_set[i]);
    
  for (i = 1; i <= Atoms; i++)
  {
    if (the_tree[i].level >= 0)
    biton(the_set[the_tree[i].level],i);
  }
}


void find_common_atom(ums_type *mol,set_type *set1[], set_type *set2[], set_type *set3)
{
  int i,j;
  set_type *temp;

  setclear(set3);
    
  temp = init_set_minbits(Atoms);

  for (i = 0;i < Atoms; i++)
  {
    for (j = 0; j <= i; j++)
    {
      setand(set1[i],set2[j],temp);
      setor(temp,set3,set3);
    }
    for (j = 0; j <= i; j++)
    {
      setand(set2[i],set1[j],temp);
      setor(temp,set3,set3);
    }
  }
  free_set(temp);

}


void 
  print_tree_set(ums_type *mol,set_type *the_set[])
{
  int i;
  char the_str[10];
  
  for (i = 0; i < Atoms; i++)
  {
    if (setcount(the_set[i]) > 0)
    {
      sprintf(the_str,"level %d",i);
      setprint(the_set[i],the_str);
    }
  }
}

int sort_rca(connect_type *a, connect_type *b)
{
  if (a->bond_order > b->bond_order)
    return 1;
  if (a->bond_order < b->bond_order)
    return -1;
  else 
    return 0;
}


void build_spanning_tree(ums_type *mol, int root, spanning_tree *tree)
{
  int i,j,k;
  

  for (i = 0; i <= Atoms; i++)
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


void build_restricted_tree(ums_type *mol, int root, int other, spanning_tree *tree)
{
  int level = 0;
  int i,j;
  int next;
  
  set_type *level_set1, *level_set2, *found_set;
  
  level_set1 = init_set_minbits(Atoms);
  level_set2 = init_set_minbits(Atoms);
  found_set = init_set_minbits(Atoms);

  tree[root].level = 0;
  tree[other].level = -1;

  biton(level_set1,root);
  biton(found_set,root);
  biton(found_set,other);

  for (i = 0; i <= Atoms; i++)
  {
    tree[i].level = 0;
    tree[i].ancestor = 0;
  }

  while (setcount(level_set1) != 0)
  {
    level++;
    setclear(level_set2);
    next = 0;
    while (next != -1)
    {
      next = NextBit(level_set1,next);
      if (next != -1)
	for (i = 0; i < Valence(next); i++)
	{
	  j = Connection(next,i);
	  if (bit_is_on(found_set,j) == FALSE)
	  {
	    tree[j].level = level;
	    tree[j].ancestor = next;
	    biton(level_set2,j);
	    biton(found_set,j);
	  }
	}
    }
    setcopy(level_set1,level_set2);
  }    
  for (i = 0; i <= Atoms; i++)
    if ((tree[i].ancestor == 0) && (i != root))
      tree[i].level = -1;
  free_set(found_set);
  free_set(level_set1);
  free_set(level_set2);
}      


void print_spanning_tree(ums_type *mol,spanning_tree *tree)
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
}

int find_closure_bonds(ums_type *mol,spanning_tree *tree, connect_type *rca)
{
  int closure;
  int i,j;
  int ring_count = 0;

  for (i = 0; i < Bonds; i++)
  {
    closure = TRUE;
    for (j = 1; j <= Atoms; j++)
    {
      if (((Start(i) == j) && (End(i) == tree[j].ancestor)) ||
	 ((End(i) == j) && (Start(i) == tree[j].ancestor)))
	closure = FALSE;
    }
    if (closure == TRUE) 
    {
/*      printf("Bond %d - %d is a ring closure bond \n",Start(i),End(i));  */
      rca[ring_count].start = Start(i);
      rca[ring_count].end = End(i);
      ring_count ++;
    }
  }
  return(ring_count);
}

void path_to_root(spanning_tree *tree, int atom, path *the_path)
{
  the_path->path_atoms[the_path->length] = atom;
  biton(the_path->path_set,atom);

  the_path->length++;
  if (tree[atom].ancestor == 0)
  {
    return;
  }
  else
  {
    path_to_root(tree,tree[atom].ancestor,the_path);
  }
}


void init_path(int atoms, path *the_path)
{
  int i;
  
  the_path->length = 0;
  the_path->bogus = FALSE;
  the_path->closure = 0;
  the_path->found = FALSE;
  for (i = 0; i < atoms; i++)
    the_path->path_atoms[i] = -1;

    setclear(the_path->path_set);
}


void print_path(path *the_path)
{
  int i;
  
  for (i = 0; i < the_path->length; i++)
    printf("%d ",the_path->path_atoms[i]);
  printf("\n");

}

void paths_to_ring(path path1, path path2, path *the_ring,int atoms)
{
  int ring_size;
  int i,j;
  
  ring_size = path1.length + path2.length - 1;
  the_ring->path_atoms = (int *)malloc(ring_size * sizeof(int));
  the_ring->length = ring_size;
  the_ring->path_set = init_set_minbits(atoms);
  j = 0;
  for (i = (path2.length - 1); i >= 0; i--)
  {
    the_ring->path_atoms[j] = path2.path_atoms[i];
    biton(the_ring->path_set,path2.path_atoms[i]);
    j++;
  }
  for (i = 1; i < path1.length; i++)
  {
    the_ring->path_atoms[j] = path1.path_atoms[i];
    biton(the_ring->path_set,path1.path_atoms[i]);
    j++;
  }
}

void show_rings(path *ring_set, int ring_count)
{
  int i,j,good_rings = 0;
  
  for (i = 0; i < ring_count; i++)
  {
    if (ring_set[i].bogus == FALSE)
    {
      for (j = 0; j < ring_set[i].length; j++)
      {
	printf("%d ", ring_set[i].path_atoms[j]);
      }
      printf("---- %d \n", ring_set[i].length);
      good_rings++;
    }
  }

  printf("Total rings = %d\n",good_rings);
}

void sort_rings(path *ring_set, int ring_count)
{
#ifndef ULTRIX
  qsort(ring_set, ring_count, sizeof(path),QSORT_PROTO comp_rings);
#else
  qsort(ring_set, ring_count, sizeof(path), comp_rings);
#endif
}

int comp_rings(path *p1, path *p2)
{
  if (p1->length > p2->length) return(1);
  if (p2->length > p1->length) return(-1);
  return(0);
}

void find_bogus_rings(path *ring_set, int ring_count,int atoms)
{
  int i,j;
  set_type *temp;
  
  temp = init_set_minbits(atoms);
    
  for (i = 0; i < ring_count; i++)
    ring_set[i].bogus = FALSE;
  for (i = 0; i < ring_count; i++)
    for (j = i; j < ring_count; j++)
    {
      if ((ring_set[i].bogus == FALSE) && (i != j))
      {
		setand(ring_set[i].path_set,ring_set[j].path_set,temp);
		if (setcount(temp) == ring_set[i].length)
		  ring_set[j].bogus  = TRUE;
      }
    }
  free_set(temp);
}  

void find_bogus_rings2(ums_type *mol,path *ring_set, int ring_count, int frj)
{
  int i,j;
  set_type *used, *same;
  int size;
  int zapped = 0;

  used = init_set_minbits(Atoms);
  same = init_set_minbits(Atoms);

  for (i = ring_count - 1; i >= 0; i--)
  {
    setclear(used);
    size = ring_set[i].length;
    for (j = 0; j < ring_count; j++)
    {
      if (ring_set[j].bogus == FALSE)
	if (ring_set[j].length <= size && i != j)
	  setor(used,ring_set[j].path_set,used);
    }
    setand(ring_set[i].path_set,used,same);

    if (ring_count - zapped == frj)
       break;

    if (setcount(same) == setcount(ring_set[i].path_set))
    {
      ring_set[i].bogus = TRUE;
      zapped++;
    }
  }

  free_set(used);
  free_set(same);
}  



void show_ring(path ring_path)
{
  int i;
  
  printf("the ring is ");
  for (i = 0; i < ring_path.length; i++)
    printf(" %d ",ring_path.path_atoms[i]);
  printf("\n");
}


int is_good_ring(path new_ring, path* ring_list, int ring_count)
{
  int i;
  set_type *temp;
  
  temp = init_set_setlen(new_ring.path_set->setlen);

  for (i = 0; i < ring_count; i++)
  {
    setand(ring_list[i].path_set,new_ring.path_set,temp);
    if (setcount(temp) == ring_list[i].length)
    {
      return(FALSE);
    }
  }
  
  free_set(temp);
  return(TRUE);
}



int build_common_array(set_type *common_set, spanning_tree *the_tree,spanning_tree *common_array)
{
  int size;
  int next = 0, k = 0;
    
  size = setcount(common_set);
  while (next != -1)
  {
    next = NextBit(common_set,next);
    if (next != -1)
    {
      common_array[k].ancestor = next;
      common_array[k].level = the_tree[next].level;
      k++;
    }
  }
  
#ifndef ULTRIX
  qsort(common_array,size,sizeof(spanning_tree),QSORT_PROTO sort_common);
#else
  qsort(common_array,size,sizeof(spanning_tree),sort_common);
#endif

  return(size);
}

int sort_common(spanning_tree *a, spanning_tree *b)
{
  if (a->level > b->level)
    return(1);
  else
    if (a->level < b->level)
      return(-1);
    else
      return(0);
}


void make_ring_ums(ums_type *mol, ums_type *new_mol, path rng)
{
  int i;
  
  new_mol->num_atoms = rng.length;
  new_mol->num_bonds = rng.length;
  initialize_ums(&new_mol);
  for (i = 1; i <= new_mol->num_atoms; i++)
  {
    strcpy(new_mol->atoms[i].type,Type(rng.path_atoms[i-1]));
    new_mol->atoms[i].point.x = X(rng.path_atoms[i-1]);
    new_mol->atoms[i].point.y = Y(rng.path_atoms[i-1]);
    new_mol->atoms[i].point.z = Z(rng.path_atoms[i-1]);
  }
  for (i = 0; i < (new_mol->num_atoms - 1); i++)
  {
    new_mol->connections[i].start = i+1;
    new_mol->connections[i].end = i+2;
  }
  new_mol->connections[new_mol->num_atoms - 1].start = new_mol->num_atoms;
  new_mol->connections[new_mol->num_atoms - 1].end = 1;
  dissect_connection_table(new_mol);
}

void
save_good_rings(ring_struct *good,ring_struct *bad,int count,int dupe_ck) 
{
  int i,j,k,*tmp;
  set_type *set;
  int duplicate;
  
  j = good->count;
  for (i = 0; i < count; i++)
  {
    duplicate = FALSE;
    for (k = 0;k < j;k++)
      if ((setcmp(good->ring_list[k].path_set,bad->ring_list[i].path_set)) && dupe_ck)
	duplicate = TRUE;

    if (!bad->ring_list[i].bogus && !duplicate)
    {
      tmp = bad->ring_list[i].path_atoms;
      set = bad->ring_list[i].path_set;
      
      bad->ring_list[i].path_atoms = NULL;
      bad->ring_list[i].path_set = NULL;
      
      good->ring_list[j] = bad->ring_list[i];
      good->ring_list[j].path_atoms = tmp;
      good->ring_list[j].path_set = set;
      j++;
    }
    else
    {
      if (bad->ring_list[i].path_atoms)
	free(bad->ring_list[i].path_atoms);
      free_set(bad->ring_list[i].path_set);
    }
  }      
  
  good->count = j;
}

