#include "bbltyp.h"

void ums_to_mini_ums(mini_ums *mini, ums_type *mol)
{
  int i;
  
  mAtoms = Atoms;
  mini->energy = Energy;

  if (!mini->atoms)
    initialize_mini(mini);
  
  for (i = 1; i <= Atoms; i++)
  {
    mX(i) = X(i);
    mY(i) = Y(i);
    mZ(i) = Z(i);
  }
}

void copy_coordinates(mini_ums *mini, ums_type *mol)
{
  int i;
  
  Energy = mini->energy;

  for (i = 1; i <= Atoms; i++)
  {
    X(i) = mX(i);
    Y(i) = mY(i);
    Z(i) = mZ(i);
  }
}


void copy_mini(mini_ums *mini1,mini_ums *mini2)
{
  int i;

  mini1->num_atoms = mini2->num_atoms;
  mini1->energy = mini2->energy;

  if (!mini1->atoms)
    initialize_mini(mini1);

  for (i = 1;i <= mini1->num_atoms;i++)
    {
      mini1->atoms[i].x = mini2->atoms[i].x;
      mini1->atoms[i].y = mini2->atoms[i].y;
      mini1->atoms[i].z = mini2->atoms[i].z;
    }
}

void initialize_mini(mini_ums *mini)
{
  int i;
  
  mini->atoms = (coord_type *)malloc((mAtoms + 1) * sizeof(coord_type));
  
  for (i = 1; i <= mAtoms; i++)
    {
      mX(i) = 0.0;
      mY(i) = 0.0;
      mZ(i) = 0.0;
    }
}

void make_mini_enantiomer(mini_ums *ena, mini_ums *mini)
{
  int i;

  ena->num_atoms = mAtoms;
  ena->energy = mEnergy;
  
  if (!ena->atoms)
    initialize_mini(ena);
  
  copy_mini(ena,mini);

  for (i = 1; i <= mAtoms; i++)
    ena->atoms[i].x *= -1;
}

void print_mini(mini_ums *mini)
{
  int i;
  
  for (i = 1; i <= mAtoms; i++)
  {
    printf("%4d %10.5f%10.5f%10.5f\n",i,mX(i),mY(i),mZ(i));
  }
}


void add_mini(mini_ums *mini1, mini_ums *mini2, mini_ums *new_mini)
{
  int i,k;
  
  new_mini->num_atoms = mini1->num_atoms + mini2->num_atoms;
  new_mini->atoms = (coord_type *)malloc((new_mini->num_atoms + 1) * sizeof(coord_type));

  for (i = 1; i <= mini1->num_atoms; i++)
  {
    new_mini->atoms[i].x = mini1->atoms[i].x;
    new_mini->atoms[i].y = mini1->atoms[i].y;
    new_mini->atoms[i].z = mini1->atoms[i].z;
  }	

  k = mini1->num_atoms;

  for (i = 1; i <= mini2->num_atoms; i++)
  {
    new_mini->atoms[i+k].x = mini2->atoms[i].x;
    new_mini->atoms[i+k].y = mini2->atoms[i].y;
    new_mini->atoms[i+k].z = mini2->atoms[i].z;
  }
}

void mini_to_set(mini_ums *mini, set_type *the_set)
{
  int i;
  
  setclear(the_set);
  for (i = 1; i <= mAtoms; i++)
  {
    if (non_zero(&mini->atoms[i]))
      biton(the_set,i);
  }
}

int non_zero(coord_type *pt)
{
  if ((pt->x == 0.0) && (pt->y == 0.0) && (pt->z == 0.0))
    return(FALSE);
  else
    return(TRUE);
}  


void adjust_mini_vector(mini_ums *mini, int core_atom, int vect_atom, double new_dist)
{
  vect_type the_vector;

  the_vector.x = mX(vect_atom) - mX(core_atom);
  the_vector.y = mY(vect_atom) - mY(core_atom);
  the_vector.z = mZ(vect_atom) - mZ(core_atom);
  
  normalize_vect(&the_vector);
  scal_x_vect(&the_vector,new_dist);
  mPoint(vect_atom) = point_plus_vector(&mPoint(core_atom),&the_vector);
}

void release_mini(mini_ums *mini)
{
  if (mini->atoms)
    free(mini->atoms);
}

void zero_mini(mini_ums *mini)
{
  int i;
  
  for (i = 1;i <= mAtoms;i++)
  {
    mX(i) = 0.0;
    mY(i) = 0.0;
    mZ(i) = 0.0;
  }
  
}

void read_mini(FILE *file,mini_ums *mini)
{
  int i;
  char buffer[100];
  
  fgets(buffer,sizeof(buffer),file);
  sscanf(buffer,"%d",&mAtoms);
  fgets(buffer,sizeof(buffer),file);
  sscanf(buffer,"%lf",&mEnergy);

  if (!mini->atoms)
  initialize_mini(mini);

  for (i = 1;i <= mAtoms;i++)
  {
    fgets(buffer,sizeof(buffer),file);
    sscanf(buffer,"%lf %lf %lf",&mX(i),&mY(i),&mZ(i));
  }
}

long int write_mini(FILE *file,mini_ums *mini)
{
  int i;

  fprintf(file,"%d\n%f\n",mAtoms,mEnergy);
  
  for (i = 1;i <= mAtoms;i++)
    fprintf(file,"%10.5f %10.5f %10.5f\n",mX(i),mY(i),mZ(i));

  return(ftell(file));
}
