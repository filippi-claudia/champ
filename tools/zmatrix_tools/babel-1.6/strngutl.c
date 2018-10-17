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
******/

#include "bbltyp.h"

/*-----------------------------------------------------------------------------
FUNCTION : rjust
PURPOSE  : right justify a string
------------------------------------------------------------------------------*/
char *rjust(char *str)
{
  int n = strlen(str);
  char *dup_str;
  
  dup_str = (char *)malloc(sizeof(char) * (strlen(str)+1));
  strcpy(dup_str,str);
  rtrim(dup_str);
  sprintf(str,"%*.*s",n,n,dup_str);
  free(dup_str);
  return(str);
}

/*-------------------------------------------------------------------------------
FUNCTION : rtrim
PURPOSE  : trim trailing spaces from a string
--------------------------------------------------------------------------------*/
char *rtrim(char *str)
{
  int n = strlen(str) - 1;
  
  while (n >= 0)
  {
    if (*(str+n) != ' ')
    {
      *(str+n+1) = '\0';
      break;
    }
    else n--;
  }
  if ((n < 0) && (str[0] == ' ')) str[0] = '\0';
  
  return str;
}

/*-------------------------------------------------------------------------------
FUNCTION : ltrim
PURPOSE  : trim leading spaces from a string
--------------------------------------------------------------------------------*/
char *ltrim(char *str)
{
  strrev(str);
  rtrim(str);
  strrev(str);
  
  return str;
}


/*---------------------------------------------------------
FUNCTION : fill_space 
PURPOSE  : pad and right justify a string
AUTHOR   : Andreas Bohne (andreas@physik.uni-hildesheim.de)
DATE     : 2-14-96
----------------------------------------------------------*/
char*
fill_space(char* dest,char* in,int gesamtl)
{
        char*   tmp     = (char*)(calloc(255,sizeof(char)));
        char*   lauf    = tmp;
        int     len     = strlen(in);
        int     ende    = gesamtl-len;
        int     i;
        for(i=0;i<ende;i++)
                *(lauf++)=' ';        /* or use '.' or '-' */
        dest = strcpy(dest,tmp);
        free(tmp);
        return dest;
}
  

/*---------------------------------------------------------
FUNCTION : strrev
PURPOSE  : reverse a string
AUTHOR   : Bob Stout
FOUND AT : ftp://www.cdrom.com/pub/algorithms/c/strrev.c
----------------------------------------------------------*/
char *strrev(char *str)
{
  char *p1, *p2;
  
      if (! str || ! *str)
	return str;
  
  for (p1 = str, p2 = str + strlen(str) - 1; p2 > p1; ++p1, --p2)
  {
    *p1 ^= *p2;
    *p2 ^= *p1;
    *p1 ^= *p2;
  }
  
  return str;
}

/*---------------------------------------------------------
FUNCTION : substr count
PURPOSE  : count how many times the string needle occurs 
           in the string haystack
DATE     : 5-21-96
----------------------------------------------------------*/

int substring_count(char *needle, char *haystack)
{
  char *pos = NULL, *newpos = NULL;
  int match_count = 0;
  int len = strlen(needle);
  
  pos = strstr(haystack,needle);
  if (pos)
  {
    match_count++;
    while (pos)
    {
      newpos = strstr(pos+(len*sizeof(char)),needle);
      if (newpos)
      {
	  match_count++;
	}
      pos = newpos;
    }
  }

  return(match_count);
}


/*--------------------------------------------------------------------
FUCTION  :     REPLACE.C
AUTHOR   :     Gilles Kohl
PURPOSE  :    Replace one string by another in a given buffer.
              This code is public domain. Use freely.
--------------------------------------------------------------------*/

char *replace_string(char *Str, char *OldStr, char *NewStr)
{
  int OldLen, NewLen;
  char *p, *q;

  if(NULL == (p = strstr(Str, OldStr)))
    return p;

  OldLen = strlen(OldStr);
  NewLen = strlen(NewStr);
  memmove(q = p+NewLen, p+OldLen, strlen(p+OldLen)+1);
  memcpy(p, NewStr, NewLen);
  return q;
}


void replace_all_string(char *Str, char *OldStr, char *NewStr)
{
  char *Start;
  
  Start = Str;
  while(NULL != (Start = replace_string(Start,OldStr,NewStr)));

}

char *ljust ( char *str )
{
  int i;
  int ilen;
  char *ptr;
      
  if (!str) return;
  if ((ilen = strlen(str)) == 0) return;

  ptr = (char *) malloc ( (ilen+1) * sizeof(char) );
  strcpy ( ptr, str );
   
  for (i = 0; i < ilen; i++) {
    if ((ptr[i] != ' ') && (ptr[i] != '\t')) break;
  }
  strcpy ( str, &ptr[i] );
  free ( ptr );

  return(str);
}            

