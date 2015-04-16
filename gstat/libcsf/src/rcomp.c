/*
 * rcomp.c
 */
#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/rcomp.c,v 2.1 1996/12/29 19:35:21 cees Exp $";
#endif


/*****************************************************************/
/*    FUNCTION :  RCOMPARE                                       */
/*****************************************************************/
/*             				                         */
/*****************************************************************/


#include "csf.h"
#include "csfimpl.h"

/* compare 2 maps for their location attributes
 * Rcompare compares 2 maps for all location attributes:
 *
 * projection,
 *
 * xUL, yUL, angle,
 *
 * cell size and
 *
 * number of rows and columns
 *
 * returns 0 if one of these attributes differ or in case of an error, 1 
 * if they are all equal.
 *
 * Merrno
 * NOT_RASTER
 */
int Rcompare(
	const MAP *m1, /* map handle 1 */
	const MAP *m2) /* map handle 2 */
{
	CHECKHANDLE_GOTO(m1, error);

	/* check if mapType is T_RASTER */
	if ((m1->main.mapType != T_RASTER)
	|| (m2->main.mapType != T_RASTER))
	{
		M_ERROR(NOT_RASTER);
		goto error;
	}

	if (
	    MgetProjection(m1) == MgetProjection(m2) && 
	    m1->raster.xUL == m2->raster.xUL &&
    	    m1->raster.yUL == m2->raster.yUL &&
	    m1->raster.cellSizeX == m2->raster.cellSizeX &&
	    m1->raster.cellSizeY == m2->raster.cellSizeY &&
	    m1->raster.angle == m2->raster.angle &&
	    m1->raster.nrRows == m2->raster.nrRows &&
	    m1->raster.nrCols == m2->raster.nrCols 
	)
		return(1);
error: 
		return(0);
}
