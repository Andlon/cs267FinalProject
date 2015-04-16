/*
 * getx0.c
$Log: getx0.c,v $
Revision 2.0  1996/05/23 13:16:26  cees
csf2clean

Revision 1.1  1996/05/23 13:11:49  cees
Initial revision

Revision 1.3  1995/11/01 17:23:03  cees
.

 * Revision 1.2  1994/09/08  17:16:23  cees
 * added c2man docs + small code changes
 *
 * Revision 1.1  1994/08/26  13:33:23  cees
 * Initial revision
 *
 */
#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/getx0.c,v 2.0 1996/05/23 13:16:26 cees Exp $";
#endif


#include "csf.h"
#include "csfimpl.h"

/* x value, upper left co-ordinate of map
 * returns the x value, upper left co-ordinate of map
 */
REAL8 RgetXUL(
	const MAP *map) /* map handle */
{
	CHECKHANDLE(map);
	return(map->raster.xUL);
}
