
#ifndef lint  
static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/strpad.c,v 2.0 1996/05/23 13:16:26 cees Exp $";
#endif

/********/
/* USES */
/********/

/* libs ext. <>, our ""  */
#include <string.h>

/* global header (opt.) and strpad's prototypes "" */
#include "csf.h"
#include "csfimpl.h"

/* headers of this app. modules called */ 

/***************/
/* EXTERNALS   */
/***************/

/**********************/ 
/* LOCAL DECLARATIONS */
/**********************/ 

/*********************/ 
/* LOCAL DEFINITIONS */
/*********************/ 

/******************/
/* IMPLEMENTATION */
/******************/
/* pad a string attribute with zeros (LIBRARY_INTERNAL)
 */
char *CsfStringPad(char *s, size_t reqSize)
{
	size_t l = strlen(s);
	PRECOND(l <= reqSize);
	(void)memset(s+l, '\0', reqSize-l);
	return s;
}
