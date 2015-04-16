/*
 * _getcell.c
 */
#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/_getcell.c,v 2.1 1996/12/29 19:35:21 cees Exp $";
#endif

#include "csf.h"
#include "csfimpl.h"

/* read one cell from a CSF raster file
 * RgetCell reads one cell value from a
 * file. 
 * returns
 * 1 if cell is successfully read,
 * 0 if not
 *
 * example
 * .so examples/csfdump1.tr
 */
size_t RgetCell(
	MAP *map,        /* map handle */
	size_t rowNr,     /* row number of cell */
	size_t colNr,     /* column number of cell */
	void *cellValue) /* write-only. buffer, large enough to hold
	                  * the value of the cell in the file and app
	                  * cell representation
	                  */
{
	return RgetSomeCells(map, 
				( (map->raster.nrCols) * rowNr) + colNr,
	             		1, cellValue);
}
