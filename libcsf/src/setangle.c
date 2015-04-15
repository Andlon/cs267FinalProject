
/*
 * setangle.c 
   $Log: setangle.c,v $
   Revision 2.0  1996/05/23 13:16:26  cees
   csf2clean

   Revision 1.1  1996/05/23 13:11:49  cees
   Initial revision

   Revision 1.2  1995/11/01 17:23:03  cees
   .

 * Revision 1.1  1994/09/07  13:23:08  cees
 * Initial revision
 *
 */
#ifndef lint  
static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/setangle.c,v 2.0 1996/05/23 13:16:26 cees Exp $";
#endif

/********/
/* USES */
/********/

/* libs ext. <>, our ""  */
#include <math.h>
#include "csf.h"
#include "csfimpl.h"

/* global header (opt.) and setangle's prototypes "" */


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

/* set the map angle cosine and sin in header (LIBRARY_INTERNAL)
 * set the map angle cosine and sin in header (LIBRARY_INTERNAL)
 * these values are only used in the co-ordinate conversion
 * routines. And since they do a counter clockwise rotation we
 * take the sine and cosine of the negative angle.
 */
void CsfSetMapAngleCosSin(
	MAP *m)		/* map handle */
{
	m->angleCos = cos(-(m->raster.angle));
	m->angleSin = sin(-(m->raster.angle));
}
