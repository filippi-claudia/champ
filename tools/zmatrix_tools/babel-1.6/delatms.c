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

FILE : delatms.c
AUTHOR(S) : Matt Stahl
DATE : 5-93
PURPOSE : remove all atoms of a specific type or types from a UMS

******/

#include "bbltyp.h"

#define NewTitle               new_mol->title
#define NewType(x)             new_mol->atoms[x].type
#define NewAtomic_number(x)    new_mol->atoms[x].atomic_number
#define NewValence(x)          new_mol->atoms[x].valence
#define NewCharge(x)           new_mol->atoms[x].charge
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

ums_type *
  delete_atoms(ums_type *mol,char *del_str)
 {  
  int heavy_count;
  char temp_str[100];
  char *pos;
  int i;

  strcpy(temp_str,del_str);
  uppercase(temp_str);
  pos = strstr(temp_str,"HOH");
  if (pos != NULL)
  {
    mol = delete_water(mol);
    for (i = 0; i < 3; i++,pos++)
      *pos = ' ';
  }
  
  if (EQ(del_str,"default"))
    strcpy(temp_str,"HC");
  
  heavy_count = tag_atoms(mol,temp_str);
  return(build_new_ums(mol,heavy_count));
}


int tag_atoms(ums_type *mol,char *del_str)
{
  int i;
  int j = 0;

  for (i = 1;i <= Atoms;i ++)
  {   
    if(strstr(del_str,Type(i)) == NULL)
    {
      j ++;
      Redo(i) = j;
    }
    else
    {
      Redo(i) = 0;
    }

  }
  return(j);
}


ums_type *
  build_new_ums(ums_type *mol,int heavy_count)
{
  int i,j = 0;
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
  strcpy(NewTitle,Title);
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
      NewCharge(Redo(i)) = Charge(i);
      NewAtomic_number(Redo(i)) = Atomic_number(i);
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
      NewBond_order(j) = Bond_order(i);
      j++;
    }
  }

  new_mol->control = mol->control;
  new_mol->next = mol->next;

  dissect_connection_table(new_mol);
  return(new_mol);
}


void dissect_connection_table(ums_type *mol)
{
  int i;
  int start, end;
  int count = 0;

  for (i = 1; i <= Atoms; i++)
    Valence(i) = 0;
  for (i = 0; i < Bonds; i++)
  {
    if (Start(i) < End(i))
    {
      start = Start(i);
      end = End(i);
    }
    else
    {
      end = Start(i);
      start = End(i);
    }
    if (!bonded(mol, start, end))
    {
      Start(count) = start;
      End(count) = end;
      Bond_order(count) = Bond_order(i);
      count++;

      Connection(start,Valence(start)) = end;
      BO(start,Valence(start)) = Bond_order(i);
      Valence(start)++;
      Connection(end,Valence(end)) = start;
      BO(end,Valence(end)) = Bond_order(i);
      Valence(end)++;
    }
  }
  Bonds = count;
}

/*
void get_inp_type(char *inp_type, ums_type *mol)
{
  int i;
  
  for (i = 0;i < MASTERSIZE;i++)
  {
    if ((master[i].type == InfileType) && (master[i].operation == input))
    {
      strcpy(inp_type,master[i].translate);
      break;
    }
  }
}
*/

void translate_del_str(char *inp_type,char *del_str)
{
  char *begin;
  char *mark;
  char std_type[40];
  char temp_str[80];
    
  if (strchr(del_str,' ') != NULL)
  {
    begin = del_str;
    
    do
    {
      mark = strchr(begin,' ');
      *mark = '\0';
      mark ++;
      get_std_type(inp_type,begin,std_type);
      strcat(temp_str,std_type);
      strcat(temp_str," ");
      begin = mark;
     
    }while (strchr(begin,' ') != NULL);
    
    get_std_type(inp_type,mark,std_type);
    strcat(temp_str,std_type);
  }
  else
  {
    get_std_type(inp_type,del_str,std_type);
    strcpy(temp_str,std_type);
  }
  strcpy(del_str,temp_str);
}






