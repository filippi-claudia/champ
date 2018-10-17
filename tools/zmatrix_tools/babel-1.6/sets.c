#include "bbltyp.h"

void
  setclear(set_type *seta)
{
  int j;
  for (j=0; j<seta->setlen; j++) 
  {
    seta->set[j] = 0;
  }
}

/*void
  biton(set_type *seta,int member)
{
   seta->set[(member / SETWORD)] |= (1 << (member % SETWORD));
}
*/

/*void
  bitoff(set_type *seta,int member)
{
  seta->set[(member / SETWORD)] &= (~(1 << (member % SETWORD)));
}
*/

void
  setcopy(set_type *seta,set_type *setb)
{
  int i;
  for (i = 0;i < seta->setlen;i++)
    seta->set[i] = setb->set[i];
}

void
  setprint(set_type *seta,char *string)
{
  int i,bit_count,line_count = 0;
  bit_count = SETWORD * seta->setlen;
  
  fprintf(stdout,"%s == [",string);
  
  for (i = 0;i < bit_count;i++)
  {
    if (bit_is_on(seta,i))
    {
      if (i < 100)
	fprintf(stdout,"%3d ",i);
      else
	fprintf(stdout,"%4d ",i);

    line_count++;
    if (line_count%15 == 0)
      fprintf(stdout,"\n");
    }
  }

  fprintf(stdout,"]\n");
}

void
  setand(set_type *seta,set_type *setb,set_type *setc)
{
  int j;

  int min = MIN(seta->setlen,setb->setlen);
  min = MIN(min,setc->setlen);

  for (j=0; j < min; j++) 
  {
    setc->set[j] = seta->set[j] & setb->set[j];
  }
}

void
  setor(set_type *seta,set_type *setb,set_type *setc)
{
  int j;
  int min = MIN(seta->setlen,setb->setlen);
  min = MIN(min,setc->setlen);
  
  for (j=0; j < min; j++) 
  {
    setc->set[j] = seta->set[j] | setb->set[j];
  }
}

void
  setorxor(set_type *seta,set_type *setb,set_type *setc)
{
  int j;
  for (j=0; j < seta->setlen; j++) 
  {
    setc->set[j] = seta->set[j] | setb->set[j];
    setc->set[j] = seta->set[j] ^ setc->set[j];
  }
}

void
  setxor(set_type *seta,set_type *setb,set_type *setc)
{
  int j;
  for (j=0; j < seta->setlen; j++) 
  {
    setc->set[j] = seta->set[j] ^ setb->set[j];
  }
}

void
  setnot(set_type *seta,set_type *setb)
{
  int j;
  for (j=0; j < seta->setlen; j++) 
  {
    seta->set[j] = ~ setb->set[j];
  }
}

int
  setcmp(set_type *seta,set_type *setb)
{
  int j;
  for (j=0; j<seta->setlen; j++) 
  {
    if((seta->set[j]) != (setb->set[j]))
    {
      return(0);
    }    
  }
  return(1);
}

int
  setissubset(set_type *seta,set_type *setb)
{
  int answer=0;
  set_type *mask;

  mask = init_set_setlen(seta->setlen);
  setand(seta,setb,mask);
  answer=setcmp(seta,mask); 
  free_set(mask);
  
  return(answer);
}

int  
  setcount(set_type *seta)
{
  int wrdcnt, shift, i=0;

    for (wrdcnt = 0; wrdcnt < seta->setlen ; wrdcnt++) 
    {
      for (shift = 0; shift < SETWORD; shift++) 
      {
	if ((seta->set[wrdcnt]>>shift)&1)
	{
	  i++;
	}
      }
    }

  return(i);
}

int 
  nextbit(set_type *seta,int last)
{
  int wrdcnt, shift;

  last++;
  
  if (last >= (seta->setlen * SETWORD))
	  return(-1);
	
  wrdcnt = last/SETWORD;
  
  for (shift = (last%32); shift < SETWORD; shift++) 
  {
      if ((seta->set[wrdcnt]>>shift)&1)
        return((SETWORD*wrdcnt) + shift);
  }
  
  for (wrdcnt = (last/SETWORD)+1; wrdcnt < seta->setlen ; wrdcnt++) 
  {
      for (shift = 0; shift < SETWORD; shift++) 
      {
	if ((seta->set[wrdcnt]>>shift)&1) 
	  return((SETWORD*wrdcnt) + shift);
      }
  }
    
  return(-1);
}

int NextBit(set_type *seta, int last)
{
   register int s,bit,wrdcnt;
   last++;

   wrdcnt = last/SETWORD;
   
   if (wrdcnt >= seta->setlen)
      return(-1);

   if (seta->set[wrdcnt] != 0)
   {
      s = (seta->set[wrdcnt]) & bitsoff[last - (wrdcnt*SETWORD)];
      if (s)
      {
         LowBit(s,bit);
         if (bit != -1)
            return(bit + (wrdcnt*SETWORD));
      }
   }
   wrdcnt++;
   
   while(wrdcnt < seta->setlen)
   {
      if (seta->set[wrdcnt] != 0)
      {
         s = seta->set[wrdcnt];
         LowBit(s, bit);

         if (bit != -1)
            return(bit+(wrdcnt*SETWORD));
      }
      wrdcnt++;
   }

   return(-1);
}

set_type *
  init_set_minbits(int minbits)
{
  set_type *new_set;
  int setlen;
  
  minbits++;
  
  setlen = (minbits/32);
  if (minbits % 32 != 0)
    setlen++;
    
  new_set = (set_type *)malloc(sizeof(set_type));
  if (new_set == NULL)
  {
    fprintf(stderr,"Unable to allocate memory for set\n");
    exit(0);
  }
  
  new_set->set = (int *)malloc(sizeof(int) * setlen);

  if (new_set->set == NULL)
  {
    fprintf(stderr,"Unable to allocate memory for set\n");
    exit(0);
  }
					     
  new_set->setlen = setlen;
  setclear(new_set);
  
  return(new_set);
}

set_type *
  init_set_setlen(int setlen)
{
  set_type *new_set;
    
  new_set = (set_type *)malloc(sizeof(set_type));

  if (new_set == NULL)
  {
    fprintf(stderr,"Unable to allocate memory for set\n");
    exit(0);
  }

  new_set->set = (int *)malloc(sizeof(int) * setlen);
  if (new_set->set == NULL)
  {
    fprintf(stderr,"Unable to allocate memory for set\n");
    exit(0);
  }
					     
  new_set->setlen = setlen;
  setclear(new_set);
  
  return(new_set);
}

set_type *
	realloc_set_setlen(set_type *set,int setlen)
{
	int i;
	if (setlen == set->setlen)
		return(set);
		
	set->set = (int *)realloc(set->set,(sizeof(int) * setlen));
	
	if (!set->set)
	{
		fprintf(stderr,"Unable to allocate memory for reallocated set\n");
		exit(0);
	}

	for (i = set->setlen;i < setlen;i++)
		set->set[i] = 0;
		
	set->setlen = setlen;
	return(set);
}

void
  free_set(set_type *set)
{
  if (set->set != NULL && set != NULL)
    free(set->set);
  
  if (set != NULL)
    free(set);
}




