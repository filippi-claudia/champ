#include "bbltyp.h"

static int first_time = TRUE;
static int mol_count = 1;

int write_csr(FILE *fp, ums_type *mol)
{
  if (first_time)
  {
    write_csr_header(fp,mol);
    first_time = FALSE;
  }
  
  write_csr_coords(fp,mol);
  mol_count++;
}

void write_csr_header(FILE *fp, ums_type *mol)
{
  char *molnames;
  int nmol, natom;

  molnames = pad_string(Title,100);

  nmol = 1;
  natom = Atoms;

  write_size(4*sizeof(char),fp);
  fwrite("V33 ",1,4*sizeof(char),fp);
  write_size(4*sizeof(char),fp);
  
  write_size(2*sizeof(int),fp);
  fwrite(&natom,1,sizeof(int),fp);
  fwrite(&nmol,1,sizeof(int),fp);
  write_size(2*sizeof(int),fp);
  
  write_size(100*sizeof(char),fp);
  fwrite(molnames,1,100*sizeof(char),fp);
  write_size(100*sizeof(char),fp);
  
  write_size(sizeof(int),fp);
  fwrite(&natom,1,sizeof(int),fp);
  write_size(sizeof(int),fp);

  free(molnames);
}



void write_csr_coords(FILE *fp, ums_type *mol)
{
  int i;
  int the_size,jconf;
  float x,y,z,energy;
  char title[100];
  char *tag;

  the_size = sizeof(int) + sizeof(float) + (80 * sizeof(char));
  
  jconf = 1;
  energy = -2.584565;

  sprintf(title,"%s:%d",Title,mol_count);
  tag = pad_string(title,80);

/*
  printf("processing  %s \n",title);
*/
  write_size(the_size,fp);
  fwrite(&jconf,1,sizeof(int),fp);
  fwrite(&energy,1,sizeof(energy),fp);
  fwrite(tag,1,80*sizeof(char),fp);
  write_size(the_size,fp);

  write_size(Atoms * sizeof(float),fp);
  for (i = 1; i <= Atoms; i++)
  {
    x = (float)X(i);
    fwrite(&x,1,sizeof(float),fp);
  }
  write_size(Atoms * sizeof(float),fp);

  write_size(Atoms *sizeof(float),fp);
  for (i = 1; i <= Atoms; i++)
  {
    y = (float)Y(i);
    fwrite(&y,1,sizeof(float),fp);
  }
  write_size(Atoms * sizeof(float),fp);

  write_size(Atoms * sizeof(float),fp);
  for (i = 1; i <= Atoms; i++)
  {
    z = (float)Z(i);
    fwrite(&z,1,sizeof(float),fp);
  }
  write_size(Atoms * sizeof(float),fp);
}

void write_size(int the_size, FILE *fp)
{
  fwrite(&the_size,1,sizeof(int),fp);
}

char *pad_string(char *input, int size)
{
  int i;
  char *output;
  
  output = (char *)malloc(size * sizeof(char));
  for (i = 0; i < size; i++)
    output[i] = ' ';
  for (i = 0; i < strlen(input); i++)
    output[i] = input[i];
  return(output);
}























