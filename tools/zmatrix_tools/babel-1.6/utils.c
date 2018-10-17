/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 

For more information please contact :

babel@mercury.aichem.arizona.edu
--------------------------------------------------------------------------------

FILE : utils.c
AUTHOR(S) : Pat Walters
DATE : 10-93
PURPOSE : General purpose utilities

******/



#include "bbltyp.h"

#ifdef MSDOS
#include <malloc.h>
#define MALLOC  _fmalloc
#define REALLOC _frealloc
#define FREE    _ffree
#else
#define MALLOC  malloc
#define REALLOC realloc
#define FREE    free
#endif

static char wstr[BUFF_SIZE];

void 
uppercase(char *str)
{
  while( *str != '\0')
  {
    if (isalpha(*str))
      if (islower(*str))
        *str = toupper(*str);
    str ++;
  }
}

void 
lowercase(char *str)
{
  while( *str != '\0')
  {
    *str = tolower(*str);
    str ++;
  }
}



void babel_init()
{
  int result;
  
  result = read_element_file();
  if (!result) 
    fatal_error("Error reading element.lis");
  read_types_table();
}

void babel_cleanup()
{
  clean_up_elements();
  clean_up_types_table();
}

void strip_return(char *the_str)
{
  int len;
  
  len = strlen(the_str) - 1;
  if (the_str[len] == '\n')
    the_str[len] = '\0';
}



int 
  print_ums(ums_type *mol)
{
  int count, i;
  char bo_char[6] = {' ','-','=','#',':',':'};
  
  fprintf(stdout,"there are %d atoms and %d bonds in this file\n",Atoms,Bonds);
  fprintf(stdout,"Title = %s\n",Title);
  
  for (count = MIN_ATOM; count <= Atoms; count++)
  {
    /*
      printf("%3datom [%c][%c][%c] %7.3f %7.3f %7.3f %5d %6.3f %5d %5d %5d %5d->",
      count,
      Type(count)[0],
      Type(count)[1],
      Type(count)[2],
      */
    printf("%3d %-6s %7.3f %7.3f %7.3f %5d %6.3f %3d %3d %3d %5d->",
	   count,
	   Type(count),
	   X(count),
	   Y(count),
	   Z(count),
	   Max_bonds(count),
	   Radius(count),
	   Redo(count),
	   count_heavy_atoms(mol, count),
	   count_free_ox(mol, count),
	   Valence(count));
    for (i = 0; i < Valence(count); i++)
    {
      printf(" %3d ",Connection(count,i));
    }
    printf("\n");
  }

  for (count = 0; count < Bonds; count++)
  {
    printf(" [%2d%c%-2d] ",Start(count),bo_char[Bond_order(count)],End(count));
    if (((count+1) % 10) == 0) printf("\n");
  }
  printf("\n");
  return(1);
}

#define MAXLINELEN 1024
      
void
show_warning ( const char *format, ... )
{
   va_list args;
   char str[MAXLINELEN];

   va_start ( args, format );
   vsprintf ( str, format, args );
   va_end ( args );

   fprintf ( stderr, "%s\n", str );
}

void
fatal_error ( const char *format, ... )
{
   va_list args;
   char str[MAXLINELEN];

   va_start ( args, format );
   vsprintf ( str, format, args );
   va_end ( args );

   fprintf ( stderr, "%s\n", str );
   exit(0);
}
#undef MAXLINELEN




int 
print_internal(ums_type *mol)
{
  int count;

  printf("there are %d atoms and %d bonds in this file\n",Atoms,Bonds);
  for (count = MIN_ATOM; count <= Atoms; count++)
  {
    printf("%4d %4s %8.3f %8.3f %8.3f %4d %4d %4d\n",
	   count,
	   Type(count),
	   mol->internal[count].r,
	   mol->internal[count].w,
	   mol->internal[count].t,
	   mol->internal[count].na,
	   mol->internal[count].nb,
	   mol->internal[count].nc);
  }
  return(1);
}
    
int 
  initialize_ums(ums_type **mol)
{
  long int byte_count;

  (**mol).num_connections_initialized = 0;
  (**mol).num_residues_initialized = 0;
  (**mol).num_internal_initialized = 0;
  
  (**mol).atoms = (atom_type HUGEPTR *)MALLOC((long int )((**mol).num_atoms + 1) * (long int)sizeof(atom_type));
  if ((**mol).atoms == NULL) 
  {
    byte_count = ((long int )((**mol).num_atoms + 1) * (long int)sizeof(atom_type));
    sprintf(wstr,"Unable to allocate memory for %d atoms requires %ld bytes\n",
	   (**mol).num_atoms,byte_count);
    if ((**mol).num_atoms > 0)
      fatal_error(wstr);
  }
  memset((**mol).atoms,0,(long int )((**mol).num_atoms + 1) * (long int)sizeof(atom_type));
  (**mol).connections = (connect_type HUGEPTR *)MALLOC((long int)(4 * (**mol).num_atoms) * (long int)sizeof(connect_type));
  if ((**mol).connections == NULL) 
  {
    byte_count = ((long int)(4 * (**mol).num_atoms) * (long int)sizeof(connect_type));
    sprintf(wstr,"Unable to allocate memory for %d bonds requires %ld bytes\n",
	    (**mol).num_atoms,byte_count);
    if ((**mol).num_bonds > 0)
      fatal_error(wstr);
  }
  memset((**mol).connections,0,(long int)((**mol).num_atoms + 1) * (long int)sizeof(connect_type));

  (**mol).num_connections_initialized = 4 * (**mol).num_atoms;

  (**mol).fract = NULL;
  (**mol).internal = NULL;
  (**mol).residues = NULL;
  (**mol).atoms[0].pos[0] = -1;
  zero_out_ums(*mol,0);
  return(1);
}

void
  zero_out_ums(ums_type *mol, int start)
{
  int i, j;
  
  Energy = 0.0;
  Next = NULL;
  
  strcpy(Title,"");
  for (i = start; i <= Atoms; i++)
  {
    X(i) = 0.0;
    Y(i) = 0.0;
    Z(i) = 0.0;
    Max_bonds(i) = 0;
    Atomic_number(i) = 0;
    for (j = 0; j < MAX_CONNECTIONS; j++)
    {
      Connection(i,j) = 0;
      BO(i,j) = 0;
    }
    Valence(i) = 0;
    BORadius(i) = 0.0;
    Radius(i) = 0.0;
    Redo(i) = 0;
    Charge(i) = 0.0;
    for (j = 0; j < 3; j++)
      mol->atoms[i].pos[j] = 0;
  }
}

    

int 
reinitialize_ums(ums_type **mol)
{
  (**mol).num_connections_initialized = 0;
  (**mol).num_residues_initialized = 0;
  (**mol).num_internal_initialized = 0;

  (**mol).atoms = (atom_type HUGEPTR *)REALLOC((**mol).atoms,((long int)(**mol).num_atoms + 1) * (long int)sizeof(atom_type));
  if ((**mol).atoms == NULL) 
  {
    printf("Unable to reallocate memory for atoms\n");
    exit(0);
  }
  
  (**mol).connections = (connect_type HUGEPTR *)REALLOC((**mol).connections,4 * (long int)(**mol).num_atoms * (long int)sizeof(connect_type));
  if ((**mol).connections == NULL) 
  {
    printf("Unable to reallocate memory for connections\n");
    exit(0);
  }
  
  if ((**mol).residues)
  {
    (**mol).residues = (res_list HUGEPTR *)REALLOC((**mol).residues,(long int)((**mol).num_atoms + 1) * (long int)sizeof(res_list));
    if ((**mol).residues == NULL) 
    {
      printf("Unable to reallocate memory for residues\n");
      exit(0);
    }
  }

  (**mol).num_connections_initialized = 4 * (**mol).num_atoms;
  (**mol).num_residues_initialized = (**mol).num_atoms + 1;
  
  return(1);
}


int 
initialize_internal(ums_type **mol)
{
  (**mol).internal = (int_type HUGEPTR *)MALLOC((long int)((**mol).num_atoms + 1) * (long int)sizeof(int_type));
  if ((**mol).internal == NULL) 
  {
    printf("Unable to allocate memory for internal coordinates\n");
    exit(0);
  }
  (**mol).num_internal_initialized = (**mol).num_atoms + 1;
  return(1);
}		

int 
initialize_fractional(ums_type **mol)
{
  (**mol).fract = (fract_type HUGEPTR *)MALLOC(sizeof(fract_type));
  if ((**mol).fract == NULL) 
  {
    printf("Unable to allocate memory for cell_parameters\n");
    exit(0);
  }
  return(1);
}		

int 
initialize_residues(ums_type **mol)
{
  (**mol).residues = (res_list HUGEPTR *)MALLOC((long int)((**mol).num_atoms + 1) * (long int)sizeof(res_list));
  if ((**mol).residues == NULL) 
  {
    printf("Unable to allocate memory for residues\n");
    exit(0);
  }
  (**mol).num_residues_initialized = (**mol).num_atoms + 1;
  return(1);
}		

int 
release_ums(ums_type *mol)
{
  (*mol).num_connections_initialized = 0;
  (*mol).num_residues_initialized = 0;
  (*mol).num_internal_initialized = 0;

  if ((*mol).atoms != NULL)
  {
    FREE((*mol).atoms);
    mol->atoms = NULL;
  }
  if ((*mol).connections != NULL)
  {
    FREE((*mol).connections);
    mol->connections = NULL;
  }
  
  if ((*mol).residues != NULL)
  {
    FREE((*mol).residues);
    mol->residues = NULL;
  }
  
  if ((*mol).internal != NULL)
  {
    FREE((*mol).internal);
    mol->internal = NULL;
  }
  
  if ((*mol).fract != NULL)
  {
    FREE((*mol).fract);
    mol->fract = NULL;
  }
  
  return(1);
}

  


void 
clean_atom_type(char id[])
{
  id[0] = toupper(id[0]);
  id[1] = tolower(id[1]);
  if (isalpha(id[1]) == 0)
    id[1] = NULL_CHAR;
  else
    id[2] = NULL_CHAR;
}   


double 
torsion(coord_type a, coord_type b, coord_type c, coord_type d)
{
  vect_type b1,b2,b3;
  vect_type c1,c2,c3;
  double l1,l2;
  vect_type v1,v2,v3,v4;
  double angle;
  
  point_to_vect(a,&v1);
  point_to_vect(b,&v2);
  point_to_vect(c,&v3);
  point_to_vect(d,&v4);

  vect_diff(&v1,&v2,&b1);
  vect_diff(&v2,&v3,&b2);
  vect_diff(&v3,&v4,&b3);
  
  cross_prod(&b1,&b2,&c1);
  cross_prod(&b2,&b3,&c2);
  cross_prod(&c1,&c2,&c3);
  
  l1 = magnitude(&c1);
  l2 = magnitude(&c2);

  if ((l1 * l2) < 0.001)
    return(0.0);
  else
  {
    angle = vect_ang(&c1,&c2);
    if (dot(&b2,&c3) > 0.0) 
      angle = angle * -1.0;
  }
  return(angle);
}

void free_line(char *the_line)
{
  if (the_line)
    free(the_line);
}

int check_for_eof(FILE *file1)
{
  int end = FALSE;
  long int position;
  char buffer[BUFF_SIZE];

  position = ftell(file1);  

  if ((fgets(buffer,sizeof(buffer),file1) == NULL))
    end = TRUE;
  else
  {
    if ((fgets(buffer,sizeof(buffer),file1) == NULL))
      if (count_tokens(buffer,"\t\n ") == 0)
	end = TRUE;
  }
  if (!end)
  {
    fseek(file1,position,0);
  }
  
  return (end);
}

void read_to_eof(FILE *file1)
{
  char buffer[BUFF_SIZE];
  
  while(fgets(buffer,sizeof(buffer),file1) != NULL);
}


  

void 
make_output_name(char *file_name, char *out_name, int counter)
{
  int indx;
  int name_length;
  int i, j = 0;
  char base[50];
  char ext[50];
  
  name_length = strlen(file_name);
  if (strchr(file_name,'.') != NULL)
  {
    indx = name_length - strlen(strchr(file_name,'.')); 
    strncpy(base,file_name,indx);
    base[indx] = '\0';
    for ( i = (indx + 1); i <= name_length; i++)
    {
      ext[j] = file_name[i];
      j ++;
    }
    sprintf(out_name,"%s%*.*d.%s",base,4,4,counter,ext);
  }
  else
  {
    strcpy(base,file_name);
    sprintf(out_name,"%s%04d",base,counter);
  }
}


  /*need to release list of mols if input doesn't match multi */

void
usage_multi_struct(void)
{
  show_warning("This appears to be a multi-structure or multi-conformer file ");
  show_warning("");
  show_warning("For multi-structure files the syntax is : ");
  show_warning("babel -i<input-type> infile all -o<output-flag> extension");
  show_warning("where all of the files will be extracted");
  show_warning("--or--");
  show_warning("babel -i<input-type> infile \"refcode ...\" -o<output-flag> extension");
  show_warning("where all of the recodes listed in the input-keywords will be extracted");
  show_warning("");
  show_warning("For multi-conformer files the syntax is : ");
  show_warning("babel -i<input-type> infile \"N1-N2\" -o<output-flag> outfile");
  show_warning("where the \"N1st\" through \"N2nd\" files will be extracted");
  show_warning("--or--");
  show_warning("babel -i<input-type> infile all -o<output-flag> outfile");
  show_warning("where all of the files will be extracted");
}


void toss(FILE *file1, int num)
{
  int i;
  char buffer[BUFF_SIZE];

  for (i = 0; i < num; i++)
    fgets(buffer,sizeof(buffer),file1);
}

int is_blank_line(char *the_line)
{
  int tokens;
  
  tokens = (count_tokens(the_line," \n\t"));
  if (tokens == 0)
    return(TRUE);
  else
    return(FALSE);
}

      
char *new_extension(char *filename,char *extension)
{
   char *last_period;

   last_period=strrchr(filename,'.');	     
   if (last_period!=NULL) *last_period=0;	
   if (extension[0]!='.') strcat(filename,".");	
   strcat(filename,extension);	  
   return(filename);
}

int is_one_three(ums_type *mol, int a, int b)
{
  int i,j;
  int conn;
  
  for (i = 0; i < Valence(a); i++)
  {
    conn = Connection(a,i);
    for (j = 0; j < Valence(conn); j++)	
      if (Connection(conn,j) == b)
	return(TRUE);
  }
  return(FALSE);
}


double out_of_plane(ums_type *mol, int atm)
{
  int c1, c2, c3, c4;  
  vect_type v,v1,sum,ck;

  c1 = atm;
  c2 = Connection(c1,0);
  c3 = Connection(c1,1);
  c4 = Connection(c1,2);

  pts_2_vect(mol,&v,c1,c2);
  normalize_vect(&v);
  pts_2_vect(mol,&v1,c1,c3);
  normalize_vect(&v1);
  vect_sum(&v,&v1,&sum);

  pts_2_vect(mol,&ck,c1,c4);

  return(180.0 - vect_ang(&sum,&ck));
}





