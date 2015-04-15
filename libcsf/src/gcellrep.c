#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/gcellrep.c,v 2.0 1996/05/23 13:16:26 cees Exp $";
#endif


#include "csf.h"
#include "csfimpl.h"

/* get cell representation
 * RgetCellRepr returns the in-file cell representation.
 * returns the cell representation 
 */
CSF_CR RgetCellRepr(
	const MAP *map) /* map handle */
{
	return(map->raster.cellRepr);
}
