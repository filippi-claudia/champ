/*****
This file is part of the Babel Program
Copyright (C) 1992-96 W. Patrick Walters and Matthew T. Stahl 
All Rights Reserved 
All Rights Reserved 
All Rights Reserved 
  
For more information please contact :
  
babel@mercury.aichem.arizona.edu
------------------------------------------------------------------------
FILE : convert.c
AUTHOR(S) : Pat Walters
DATE : 4-93
PURPOSE : main program
******/

#include "bbltyp.h"
#include "bblmast.h"

#undef MAC 

#ifdef MAC
#include <console.h>
#endif


static warning wstr;

static char *program_name;
int use_title = FALSE;

#ifdef MSDOS
#define UPPER_LIMIT 10000
#else
#define UPPER_LIMIT 1000000
#endif

#ifndef LIBRARY
int main(int argc, char *argv[])
{
  ums_type *mol = NULL;
  int need_to_do_outputs = FALSE;

#ifdef MAC
  argc = ccommand(&argv);
#endif

  babel_init();
  program_name = argv[0];
  if (argc == 1)
  {
    usage();
    exit(0);
  }    
  
  mol = (ums_type *)malloc(sizeof(ums_type));
  mol->control = (bbl_control *)malloc(sizeof(bbl_control));
  if (mol == NULL)
    fatal_error("Could not allocate memory for UMS");
  
  init_babel_control(mol);
  Size = MASTERSIZE;
  
  process_command_line_args(argc,argv,mol);

  if ((InfileType == none) && (UseMenus == FALSE))
    fatal_error("No input flag specified");

  if ((OutfileType == none) && (UseMenus == FALSE))
    fatal_error("No output flag specified");
  
  if (Spline)
  {
    do_spline(mol);
    return(0);
  }
  
  if (UseMenus)
    babel_menus(mol);
  else
  {
    mol = do_inputs(mol);
  }

  if (mol->control)
    free(mol->control);
  if (mol)
    free(mol);

  return(0);
}
#endif

void new_control(ums_type *mol)
{

  mol->control = (bbl_control *)malloc(sizeof(bbl_control));

  if (!mol->control)
    fatal_error("Unable to allocate memory for babel control");
  
  init_babel_control(mol);
}


void init_babel_control(ums_type *mol)
{
  strcpy(InfileName,"");
  strcpy(OutfileName,"");
  strcpy(InputKeywords,NOKEY);
  strcpy(OutputKeywords,NOKEY);    
  strcpy(DeleteStr,"");         
  InfileType = none;
  OutfileType = none;
  UseMenus = FALSE;
  Verbose = FALSE;
  AddHydrogens = FALSE;
  DeleteAtoms = FALSE;
  ReaderFunction = NULL;
  WriterFunction = NULL;
  LowerLimit = 1;
  UpperLimit = UPPER_LIMIT;
  Multi = single_struct;
  Spline = FALSE;
  Align = FALSE;
  CenterMol = FALSE;
  Increment = 0;
  DoRenum = FALSE;
  Precipitate = FALSE;
  MakeNewFile = TRUE;
  CalcCharges = FALSE;  
  StdOrientation = FALSE;
  PushHydrogens = FALSE;
  NoDummy = FALSE;
}

void 
  process_command_line_args(int argc, char *argv[],ums_type *mol)
{
  char tempstr[BUFF_SIZE], target[BUFF_SIZE];
  
  while (argc > 1) 
  {
    if (argv[1][0] == '-')
    {
      switch (argv[1][1])
      {
      case 'a' :
	lowercase(&argv[1][1]);
	if EQ(&argv[1][1],"align")
	  Align = TRUE;
	break;

      case 'b' :
	if (EQ(&argv[1][1],"blurb"))
	  write_blurb();
	if (EQ(&argv[1][1],"bbldef"))
	  write_bbldef();
	exit(0);
	break;

      case 'c' :
	lowercase(&argv[1][1]);
	if (EQ(&argv[1][1],"center"))
	  CenterMol = TRUE;
	if (EQ(&argv[1][1],"charge"))
	  CalcCharges = TRUE;
	break;

      case 'd':
	DeleteAtoms = TRUE;
	if ((argv[2] != NULL) && (argv[2][0] != '-'))
	  strcpy(DeleteStr,argv[2]);
	else 
	{
	  sprintf(wstr,"Default type to be deleted is hydrogen");
	  show_warning(wstr);
	  strcpy(DeleteStr,"default");
	}
	break;

      case 'h':
	lowercase(&argv[1][1]);
	if (EQ(&argv[1][1],"help"))
	{
	  if (argv[2] == NULL)
	    usage();
	  else
	  {
	    if (EQ(argv[2],"inputs"))
	      show_inputs();
	    else
	      if (EQ(argv[2],"outputs"))
		show_outputs();
	  }
	  exit(0);
	}
	else
	AddHydrogens = TRUE;
	break;

      case 'i':
	translate_input_code(&argv[1][2],mol);
	if (InfileType == none)
	{
	  sprintf(wstr,"%s is not a supported input file type \n",argv[1]);
	  fatal_error(wstr);
	}

	if (argv[2][0] == '-')
	  strcpy(InfileName,"CON");
	else
	  strcpy(InfileName,argv[2]);

	if ((argv[3] != NULL) && (argv[3][0] != '-'))
	  strcpy(InputKeywords,argv[3]);
	else 
	  strcpy(InputKeywords,"KEYWORDS GO HERE");
	set_limits(mol);
	break;

      case 'm':
	Verbose = TRUE;
	UseMenus = TRUE;
	break;

      case 'n' :
	lowercase(&argv[1][1]);
	if (EQ(&argv[1][1],"nodummy"))
	  NoDummy = TRUE;
	break;

      case 'o':
	if (EQ(&argv[1][1],"orient"))
	{
	  StdOrientation = TRUE;
	}
	else
	{
	  translate_output_code(&argv[1][2],mol);
	  if (OutfileType == none)
	  {
	    sprintf(wstr,"%s is not a supported output file type \n",argv[1]);
	    fatal_error(wstr);
	  }
	  if (argv[2] != NULL)
	  {
	    strcpy(OutfileName,argv[2]);
	    strcpy(BaseName,OutfileName);
	    strcpy(tempstr,OutfileName);
	    uppercase(tempstr);
	    

	    get_token(target,tempstr,". \n\t",1);
	    
	    if (EQ(target,"TITLE"))
	    {
	      Multi = multi_struct;
	      use_title = TRUE;
	      get_token(target,OutfileName,". \n\t",2);
	      if (target != NULL)
		strcpy(DefaultExtension,target);
	      else
		DefaultExtension[0] = '\0';
	    }
	    else
	      Multi = multi_conf;
	  }
	  else
	    strcpy(OutfileName,"CON");
	  if ((argv[3] != NULL) && (argv[3][0] != '-'))
	    strcpy(OutputKeywords,argv[3]);
	  else 
	    strcpy(OutputKeywords,"KEYWORDS GO HERE");
	}
	break; /* case 'o' */

      case 'p' :
	if (EQ(&argv[1][1],"precip"))
	  Precipitate = TRUE;
	if (EQ(&argv[1][1],"push"))
	  PushHydrogens = TRUE;
	break;

      case 'r' :
	lowercase(&argv[1][1]);
	if EQ(&argv[1][1],"renum")
	  DoRenum = TRUE;
	if ((argv[2] != NULL) && (argv[2][0] != '-'))
	{
	  if (isdigit(argv[2][0]))
	    NewBase = atoi(argv[2]);
	  else
	    NewBase = 1;
	}
	else
	  NewBase = 1;
	break;

      case 's':
	if (EQ(&argv[1][1],"spline"))
	{
	  Spline = TRUE;
	  if ((argv[2] != NULL) && (argv[2][0] != '-'))
	    Increment = atoi(argv[2]) + 2;
	}
	if (EQ(&argv[1][1],"split"))
	{
	  Multi = multi_struct;
	}
	break;

      case 'v':
	lowercase(&argv[1][1]);
	if (EQ(&argv[1][1],"version"))
	{
	  printf("Babel version %s %s %s \n",BABEL_VERSION,__DATE__,__TIME__);
	  exit(0);
	}
	else
	  Verbose = TRUE;
	break;
      }
    }
    argv++;
    argc--;
  }
}

void translate_input_code(char *code, ums_type *mol)
{
  int i;
  
  for (i = 0; i < MASTERSIZE; i++)
    if (EQ(code,master[i].code) && (master[i].operation == input))
    {
      InputInfo = master[i];
    }
}

void translate_output_code(char *code,ums_type *mol)
{
  int i;
  for (i = 0; i < MASTERSIZE; i++)
    if (EQ(code,master[i].code) && (master[i].operation == output))
    {
      OutputInfo = master[i];
    }
}

ums_type *do_inputs(ums_type *mol)
{		
  FILE *infile, *outfile;
  int end = FALSE, result = 0;
  int in_count = 1;
  int out_count = 1;

  if (EQ(InfileName,OutfileName) && NOTEQ(InfileName,"CON"))
  {
    fatal_error("FATAL ERROR : Input and Output file names are the same !!");
  }
  
  if (Verbose)
  {
    sprintf(wstr,"Reading %s file %s",InputTypeName,InfileName);
    show_warning(wstr);
  }
  
  if (NOTEQ(InfileName,"CON"))
    infile = open_read(InfileName);
  else
    infile = stdin;
  
  if ((Multi != multi_struct) && (UseMenus == FALSE))
    outfile = open_write(OutfileName);
  
  while (!end)
  {
    if (in_count > UpperLimit)
    {
      end = TRUE;
      break;
    }

    result = ReaderFunction(infile,mol);
    
    end = check_for_eof(infile);    
    
    if (UseMenus)
    {
      fclose(infile);
      return(mol);
    }

    if (want_this_file(mol,in_count,end))
    {
      if (Verbose)
      {
	fprintf(stderr,"Structure %5d ",in_count,Multi);
      }


      if (Multi == multi_struct)
      {
	puts("Multi == multi_struct");
	generate_outfile_name(mol,out_count);
	outfile = open_write(OutfileName);
	do_outputs(outfile,mol);
	close_file(OutfileName,outfile);
      }
      else
      {
	mol = do_outputs(outfile,mol);
      }

      out_count++;
    }
    in_count++;
  }
  
  if (Multi != multi_struct)
    close_file(OutfileName,outfile);
  if (infile)
    close_file(InfileName,infile);

  return(mol);
}

ums_type *do_outputs(FILE *outfile, ums_type *mol)
{		
  int result = 0;
  vect_type v;
  ums_type *temp;

  
  if (AddHydrogens)
  {
    if (OutfileType == smiles)
      add_2d_hydrogens(mol);
    else
      add_hydrogens(mol);
  }

  if (CalcCharges)
  {
    calc_gasteiger_charges(mol);
  }
  
  if (DoRenum)
  {
    temp = renum_for_zmat(mol,NewBase);
    temp->control = mol->control;
    release_ums(mol);
    free(mol);
    mol = temp;
  }
  
  if (CenterMol)
    center_at_origin(mol,&v);

  if (Align)
    AlignMol(mol);

  if (StdOrientation)
    set_std_orientation(mol);
  
  if (Atoms > 0)
  {  
    if (DeleteAtoms)
    {
      temp = delete_atoms(mol,DeleteStr);
      temp->control = mol->control;
      release_ums(mol);
      free(mol);
      mol = temp;
    }
    
    if (Precipitate)
    {
      temp = precipitate(mol);
      temp->control = mol->control;
      release_ums(mol);
      free(mol);
      mol = temp;
    }

    if (PushHydrogens)
    {
      push_hydrogens_to_end(mol);
      temp = build_new_ums(mol,Atoms);
      temp->control = mol->control;
      release_ums(mol);
      free(mol);
      mol = temp;
    }

    if (Title[0] == '\0')
      strcpy(Title,OutfileName);

    if (Verbose)
    {
      fprintf(stderr,"%s Writing %d Atoms and %d Bonds to %s file %s\n",
	      Title,Atoms,Bonds,OutfileTypeName,OutfileName); 
    }
    if (!(NoDummy) || dummy_check(mol))
    {
      result = WriterFunction(outfile,mol);
    }
    else
    {
      fprintf(stderr,"%s rejected because of dummy atoms\n",Title);
    }
  }
  else
  {
    show_warning("No atoms found in this structure");
  }

  release_ums(mol);
  return(mol);
}

void generate_outfile_name(ums_type *mol,int count)
{
  if (NOTEQ(OutfileName,"CON"))
  {
    if (Multi == multi_struct)
    {
      if (use_title)
      {
	strcpy(Title,trim_spaces(Title));
	if (DefaultExtension[0] != '\0')
	  sprintf(OutfileName,"%s.%s",Title,DefaultExtension);
	else
	  sprintf(OutfileName,"%s",Title);
      }
      else
	if ((UpperLimit - LowerLimit) > 0)
	  make_output_name(BaseName,OutfileName,count);
    }
    else
      if ((UpperLimit - LowerLimit) > 0)
      	make_output_name(BaseName,OutfileName,count);
    lowercase(OutfileName);
  }
}

void set_limits(ums_type *mol)
{
  char token1[25],token2[25];
  
  uppercase(InputKeywords);
  if (EQ(InputKeywords,"ALL"))
  {
    LowerLimit = 1;
    UpperLimit = UPPER_LIMIT - 1;
  }
  else
    if (strchr(InputKeywords,'-'))
    {
      get_token(token1,InputKeywords,"-",1);
      get_token(token2,InputKeywords,"-",2);
      
      if ((isdigit(token1[0])) && (isdigit(token2[0])))
      {
	LowerLimit = atoi(token1);
	UpperLimit = atoi(token2);
      }
    }
}


int want_this_file(ums_type *mol, int counter, int end)
{

  /* single structure -- write it out */
  if ((counter == 1) && (end == 1))
  {
    Multi = single_struct;
    return(TRUE);
  }

  /* no keywords -- do all structures */
  if ((EQ(InputKeywords,NOKEY)) || (EQ(InputKeywords,"ALL")))
    return(TRUE);

  /* just do the last structure */
  if (EQ(InputKeywords,"LAST"))
    if (end)
    {
      return(TRUE);
    }
    else
      return(FALSE);

  if ((strlen(Title) > 0) && (strstr(InputKeywords,Title) != NULL))
    return(TRUE);

  switch(Multi)
  {
  case multi_conf :
    if ((counter >= LowerLimit) && (counter <= UpperLimit))
      return(TRUE);
    break;
  case multi_struct :
    return(TRUE);
    break;
  default :
    return(FALSE);
  }

  return(FALSE);
}

void 
  usage()
{
  sprintf(wstr,"Babel %s -- %s -- %s",BABEL_VERSION,__DATE__,__TIME__);
  puts(wstr);
  sprintf(wstr,"for menus type -- %s -m\n",program_name);
  puts(wstr);
  sprintf(wstr,"Usage is : \n%s [-v] -i<input-type> <name> -o<output-type> <name> \"<keywords>\"\n",program_name);
  puts(wstr);
  sprintf(wstr,"\nCurrently supported input types\n");
  puts(wstr);
  show_inputs();
  sprintf(wstr,"\nCurrently supported output types\n");
  puts(wstr);
  show_outputs();
}  

void write_bbldef()
{
  int i;
  char bbldef[256];
  FILE *fp;

    if (getenv("HOME"))
      sprintf(bbldef,"%s/.bbldef",getenv("HOME"));
    else
      fatal_error("Please define HOME environment variable");

  fp = open_write(bbldef);
  for (i = 0; i < MASTERSIZE; i++)
  {
    if (master[i].operation == input)
      fprintf(fp,"input %s\n",master[i].type_name);
    if (master[i].operation == output)
      fprintf(fp,"output %s\n",master[i].type_name);
  }

  fprintf(fp,"default input MDL Isis\n");
  fprintf(fp,"default output Sybyl Mol2\n");

  if (fp)
    fclose(fp);
}





void show_inputs()
{
  int i;
  
  for (i = 0; i < MASTERSIZE; i++)
    if (master[i].operation == input)
    {
      sprintf(wstr,"\t%s -- %s file",master[i].code,master[i].type_name);
      puts(wstr);
    }
}

void show_outputs()
{
  int i;

  for (i = 0; i < MASTERSIZE; i++)
    if (master[i].operation == output)
    {
      sprintf(wstr,"\t%s -- %s file",master[i].code,master[i].type_name);
      puts(wstr);
    }
}

int output_all_formats(FILE *fp,ums_type *mol)
{
  int i;
  int result;
   
  for (i = 0; i < MASTERSIZE; i++)
    if ((master[i].operation == output) && (master[i].type != diagnostics) && (master[i].type != gaussian_template))
    {
      fprintf(stdout,"\n%s FILE\n",master[i].type_name); 
      strcpy(Title,"TITLE");
      result = master[i].func(stdout,mol);
    }
  return(TRUE);
}


void write_blurb()
{
  int i,j;
  
  printf("Babel will read the following file types :\n");  
  j = 0;
  for (i = 0; i < MASTERSIZE; i++)
    if (master[i].operation == input)
    {
      j++;
      printf("%-25s",master[i].type_name);
      if (((j % 3) == 0) && j > 0)
	printf("\n");
    }
  printf("\n");
  
  printf("\nBabel will write the following file types :\n");
  j = 0;
  for (i = 0; i < MASTERSIZE; i++)
    if ((master[i].operation == output) && (master[i].type != diagnostics))  
    {
      j++;
      printf("%-25s",master[i].type_name);
      if (((j % 3) == 0) && j > 0)
	printf("\n");
    }
  printf("\n");
}



