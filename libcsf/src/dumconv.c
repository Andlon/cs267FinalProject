
/*
 * dumconv.c 
   $Log: dumconv.c,v $
   Revision 2.0  1996/05/23 13:16:26  cees
   csf2clean

   Revision 1.1  1996/05/23 13:11:49  cees
   Initial revision

   Revision 1.2  1995/11/01 17:23:03  cees
   .

 * Revision 1.1  1994/09/09  12:17:59  cees
 * Initial revision
 *
 */
#ifndef lint  
static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/dumconv.c,v 2.0 1996/05/23 13:16:26 cees Exp $";
#endif

/********/
/* USES */
/********/

/* libs ext. <>, our ""  */
#include "csf.h"
#include "csfimpl.h"

/* global header (opt.) and dumconv's prototypes "" */

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

/* ARGSUSED */

/* dummy conversion (LIBRARY_INTERNAL)
 * does nothing
 */
void CsfDummyConversion(
	size_t nrCells,  
	void  *buf)
{
	/* nothing */
}
