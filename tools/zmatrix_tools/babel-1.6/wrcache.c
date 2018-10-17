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

FILE : wrcache.c
AUTHOR(S) : Jack Smith (Union Carbide)
DATE : 02-09-94
PURPOSE : Routines to write a CAChe MolStruct file
******/

#include "bbltyp.h"

int 
write_cache(FILE *file1, ums_type *mol)
{ 
  int i,k,l;
  char type_name[5];
  int result;
  int atomic_number;
  int formal_charge;
  static char bondType[4][7] = {"single","double","triple","weak"};

  
  assign_bond_order(mol);

  fprintf(file1,"molstruct88_Apr_30_1993_11:02:29 <molecule> 0x1d00\n");
  fprintf(file1,"Written by Molecular Editor on <date>\n");
  fprintf(file1,"Using data dictionary         9/9/93  4:47 AM\n");
  fprintf(file1,"Version 6\n");
  fprintf(file1,"local_transform\n");
  fprintf(file1,"0.100000 0.000000 0.000000 0.000000\n");
  fprintf(file1,"0.000000 0.100000 0.000000 0.000000\n");
  fprintf(file1,"0.000000 0.000000 0.100000 0.000000\n");
  fprintf(file1,"0.000000 0.000000 0.000000 1.000000\n");
  fprintf(file1,"object_class atom\n");
  fprintf(file1,"property xyz_coordinates MoleculeEditor angstrom 6 3 FLOAT\n");
  fprintf(file1,"property anum MoleculeEditor unit 0 1 INTEGER\n");
  fprintf(file1,"property sym MoleculeEditor noUnit 0 2 STRING\n");
  fprintf(file1,"property chrg MoleculeEditor charge_au 0 1 INTEGER\n");
  fprintf(file1,"property rflag MoleculeEditor noUnit 0 1 HEX\n");
  fprintf(file1,"ID xyz_coordinates             anum sym	chrg rflag\n");

  for(i = 1;i <= Atoms; i++)
  {
    result = get_output_type(i,"XYZ",Type(i),type_name,all_caps);
    formal_charge = 0;
    atomic_number = SymbToNum(type_name);
    
    fprintf(file1,"%3d %10.6f %10.6f %10.6f %2d %2s %2d 0x7052\n",
            i,X(i),Y(i),Z(i),atomic_number,type_name,formal_charge);
  }

  fprintf(file1,"property_flags:\n");
  fprintf(file1,"object_class bond\n");
  fprintf(file1,"property rflag MoleculeEditor noUnit 0 1 HEX\n");
  fprintf(file1,"property type MoleculeEditor noUnit 0 1 NAME\n");
  fprintf(file1,"property bond_order MoleculeEditor noUnit 4 1 FLOAT\n");
  fprintf(file1,"ID rflag type bond_order\n");

  for(i = 0; i < Bonds; i++)
  {
    fprintf(file1,"%3d 0x7005 %s\n",
            i+1,bondType[Bond_order(i)-1]);
  }

  fprintf(file1,"property_flags:\n");
  fprintf(file1,"object_class connector\n");
  fprintf(file1,"property dflag MoleculeEditor noUnit 0 1 HEX\n");
  fprintf(file1,"property objCls1 MoleculeEditor noUnit 0 1 NAME\n");
  fprintf(file1,"property objCls2 MoleculeEditor noUnit 0 1 NAME\n");
  fprintf(file1,"property objID1 MoleculeEditor noUnit 0 1 INTEGER\n");
  fprintf(file1,"property objID2 MoleculeEditor noUnit 0 1 INTEGER\n");
  fprintf(file1,"ID dflag objCls1 objCls2 objID1 objID2\n");

  for (k=0,l=1; k<Bonds; k++)                          
  {
    fprintf(file1,"%3d 0xa1 atom bond %d %d\n",l++,Start(k),k+1);
    fprintf(file1,"%3d 0xa1 atom bond %d %d\n",l++,End(k),k+1);
  }

  fprintf(file1,"property_flags:\n");
  return(TRUE);
}

int	SymbToNum(char	*atSymb)
{
	char	char1,char2;

	char1 = *atSymb;
	if (islower(char1)) char1 = toupper(char1);
	char2 = *(atSymb+1);
	if (isupper(char2)) char2 = tolower(char2);

	switch(char1)
	{
		case 'A':
		switch(char2)
		{
			case 'c': return(89); /* Actinium */
			case 'g': return(47); /* Silver */
			case 'l': return(13); /* Aluminum */
			case 'r': return(18); /* Argon */
			case 's': return(33); /* Arsenic */
			case 't': return(85); /* Astatine */
			case 'u': return(79); /* Gold */
		}

		case 'B':
		switch(char2)
		{
			case 'a': return(56); /* Barium */
			case 'e': return(4);  /* Beryllium */
			case 'i': return(83); /* Bismuth */ 
			case 'r': return(35); /* Bromine */
			default:  return(5);  /* Boron */
		}

		case 'C':
		switch(char2)
		{
			case 'a': return(20); /* Calcium */
			case 'd': return(48); /* Cadmium */
			case 'e': return(58); /* Cerium */
			case 'f': return(98); /* Californium */
			case 'l': return(17); /* Chlorine */
			case 'o': return(27); /* Cobalt */
			case 'r': return(24); /* Chromium */
			case 's': return(55); /* Cesium */
			case 'u': return(29); /* Copper */
			default:  return(6);  /* Carbon */
		}

		case 'D':
		switch(char2)
		{
			case 'y': return(66); /* Dysprosium */
			default:  return(1);  /* Deuterium */
		}

		case 'E': break;

		case 'F':
		switch(char2)
		{
			case 'e': return(26); /* Iron */
			case 'r': return(87); /* Francium */
			default:  return(9);  /* Fluorine */
		}

		case 'G':
		switch(char2)
		{
			case 'a': return(31); /* Gallium */
			case 'e': return(32); /* Germanium */
		}

		case 'H':
		switch(char2)
		{
			case 'e': return(2);  /* Helium */
			case 'f': return(72); /* Hafnium */
			case 'g': return(80); /* Mercury */
			default:  return(1);  /* Hydrogen */
		}

		case 'I':
		switch(char2)
		{
			case 'e': return(2);  /* Helium */
			case 'n': return(49); /* Indium */
			case 'r': return(77); /* Iridium */
			default:  return(53); /* Iodine */
		}

		case 'J': break;

		case 'K':
		switch(char2)
		{
			case 'r': return(36); /* Krypton */
			default:  return(19); /* Potassium */
		}

		case 'L':
		switch(char2)
		{
			case 'a': return(57); /* Lanthanum */
			case 'i': return(3);  /* Lithium */
		}

		case 'M':
		switch(char2)
		{
			case 'g': return(12); /* Magnesium */
			case 'n': return(25); /* Manganese */
			case 'o': return(42); /* Molybdenum */
		}

		case 'N':
		switch(char2)
		{
			case 'a': return(11); /* Sodium */
			case 'b': return(41); /* Niobium */
			case 'e': return(10); /* Neon */
			case 'i': return(28); /* Nickel */
			default:  return(7);  /* Nitrogen */
		}

		case 'O':
		switch(char2)
		{
			case 's': return(76); /* Osmium */
			default:  return(8);  /* Oxygen */
		}

		case 'P':
		switch(char2)
		{
			case 'b': return(82); /* Lead */
			case 'd': return(46); /* Palladium */
			case 'o': return(84); /* Polonium */
			case 't': return(78); /* Platinum */
			default:  return(15); /* Phosphorus */
		}

		case 'Q': break;

		case 'R':
		switch(char2)
		{
			case 'a': return(88); /* Radium */
			case 'b': return(37); /* Rubidium */
			case 'e': return(75); /* Rhenium */
			case 'h': return(45); /* Rhodium */
			case 'n': return(86); /* Radon */
			case 'u': return(44); /* Ruthenium */
		}

		case 'S':
		switch(char2)
		{
			case 'b': return(51); /* Antimony */
			case 'c': return(21); /* Scandium */
			case 'e': return(34); /* Selenium */
			case 'i': return(14); /* Silicon */
			case 'n': return(50); /* Tin */
			case 'r': return(38); /* Strontium */
			default:  return(16); /* Sulfur */
		}

		case 'T':
		switch(char2)
		{
			case 'c': return(43); /* Technetium */
			case 'e': return(52); /* Tellurium */
			case 'i': return(22); /* Titanium */
			case 'h': return(90); /* Thorium */
			case 'l': return(81); /* Thallium */
		}

		case 'U':         return(92); /* Uranium */
		case 'V':         return(23); /* Vanadium */
		case 'W':         return(74); /* Wolfram (Tungsten) */
		case 'X':         return(54); /* Xenon */
		case 'Y':         return(39); /* Yttrium */

		case 'Z':
		switch(char2)
		{
			case 'n': return(30); /* Zinc */
			case 'r': return(40); /* Zirconium */
		}
	}

	return(0);
}
