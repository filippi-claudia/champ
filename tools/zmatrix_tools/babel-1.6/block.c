/*
   block.c

   allows multiple arrays to be malloc'ed in one easy step

   Simon Kilvington, University of Southampton, 1995
*/

#include "bbltyp.h"

/* the routine that does all the work */
static int block__doalloc(int, block_ptr *, const char *, va_list);

/*
   block_alloc
   allocates a block of memory for a set of arrays
   the types string determines what the arrays are arrays of, one char per array, types are:
   c - char
   i - int
   f - float
   d - double
   B - int
   v - void *
   eg int *intarray;
      float *floatarray1, *floatarray2;
      block_alloc(&handle, "iff", &intarray, nentries1, &floatarray1, nentries2, &floatarray2, nentries3);
      ...code...
      block_free(&handle);
   returns TRUE if all went well
*/

int
block_alloc(block_ptr *handle, const char *types, ...)
{
   int okay;
   va_list ap;

   va_start(ap, types);
   okay = block__doalloc(FALSE, handle, types, ap);
   va_end(ap);

   return okay;
}

/*
   block_calloc
   as block_alloc, but all the space is initialised to 0
*/

int
block_calloc(block_ptr *handle, const char *types, ...)
{
   int okay;
   va_list ap;

   va_start(ap, types);
   okay = block__doalloc(TRUE, handle, types, ap);
   va_end(ap);

   return okay;
}

/*
   block_free
   deallocates space claimed with block_[c]alloc
*/

void
block_free(block_ptr *handle)
{
   free(*handle);
   *handle = NULL;

   return;
}

/*
   block__doalloc
   does all the work for both the above alloc'ing routines
*/

#define typesize(CHAR, TYPE)	if(types[i] == CHAR)	{			\
				   (void) va_arg(ap, TYPE **);			\
				   size += va_arg(ap, int) * sizeof(TYPE);	\
				   continue;		}

#define typeptr(CHAR, TYPE)	if(types[i] == CHAR)	{					\
				   array = (void *) va_arg(ap, TYPE **);			\
				   *((TYPE **) array) = (TYPE *) ((long) (*handle) + ptr);	\
				   ptr += va_arg(ap, int) * sizeof(TYPE);			\
				   continue;		}

static int
block__doalloc(int clear, block_ptr *handle, const char *types, va_list initap)
{
   va_list ap;
   int i, size;
   long ptr;
   void *array;

/* calc how much space we are gonna need */
   ap = initap;
   size = 0;
   for(i=0; types[i] != '\0'; i++)
   {
      typesize('c', char);
      typesize('i', int);
      typesize('f', float);
      typesize('d', double);
      typesize('B', int);
      typesize('v', void *);
   }

   *handle = (block_ptr) ((clear) ? calloc(size, 1) : malloc(size));

/* set up the ptrs if we can alloc the memory */
   if(*handle != NULL)
   {
      ap = initap;
      ptr = 0;
      for(i=0; types[i] != '\0'; i++)
      {
	 typeptr('c', char);
	 typeptr('i', int);
	 typeptr('f', float);
	 typeptr('d', double);
	 typeptr('B', int);
	 typeptr('v', void *);
      }
   }

   return (*handle != NULL);
}

