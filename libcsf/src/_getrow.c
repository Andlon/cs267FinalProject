/*
 * _getrow.c
 */
#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/_getrow.c,v 2.0 1996/05/23 13:16:26 cees Exp $";
#endif

#include "csf.h"
#include "csfimpl.h"

/* read one row from a CSF raster file
 * RgetRow reads one row of cells from a
 * file. 
 * returns
 * Number of cells successfully read
 *
 * example
 * .so examples/_row.tr
 */
size_t RgetRow(
	MAP *map,        /* map handle */
	size_t rowNr,     /* row number to be read */
	void *buf)       /* write-only. buffer large enough to hold
	                  * cell values of one row in both the file
	                  * and in-app cell representation
	                  */
{
	return RgetSomeCells(map,
				map->raster.nrCols*rowNr,
				map->raster.nrCols, 
				buf) ;
}
