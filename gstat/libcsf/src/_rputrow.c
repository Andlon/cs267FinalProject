/*
 * _rputrow.c
 */
#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/_rputrow.c,v 2.0 1996/05/23 13:16:26 cees Exp $";
#endif


#include "csf.h"
#include "csfimpl.h"


/* write one row to a CSF raster file
 * RputRow writes one row of cell values to a
 * file.
 * returns
 * number of cells successfully written. Should be equal
 * to the number of columns if everything is OK.
 *
 * example
 * .so examples/_row.tr
 */
size_t RputRow(
MAP *map,         /* map handle */
size_t rowNr,      /* Row number of row */
void *buf)        /* read-write. Buffer large enough to
                   * hold one row in the in-file cell representation
                   * or the in-app cell representation.
                   * If these types are not equal then the buffer is
                   * converted from the in-app to the in-file 
                   * cell representation. 
                   */
{
	return RputSomeCells(map, (map->raster.nrCols)*rowNr,
 	                             map->raster.nrCols, buf) ;
}
