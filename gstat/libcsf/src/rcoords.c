/*
 * rcoords.c
 */
#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/rcoords.c,v 2.1 1996/12/29 19:35:21 cees Exp $";
#endif

#include "csf.h"
#include "csfimpl.h"

/* compute true world co-ordinate of a pixel
 * RrowCol2Coords computes the true world co-ordinate from a 
 * row, column index. 
 * The row, column co-ordinate 
 * don't have to be on the map. They are just relative to upper left position.
 * For example (row,col) = (-1,0) computes the (x,y) co-ordinate of
 * the pixel that is right above upper left pixel. 
 *
 * returns
 *  0  if the co-ordinate is outside the map.
 *  1 if inside.
 * -1 in case of an error.
 *
 * Merrno
 * ILL_CELLSIZE
 */
int RgetCoords(
	const MAP *m,	/* map handle */
	int inCellPos,  /* nonzero if you want the co-ordinate
			 * at the centre of the cell, 0 if you
			 * want the upper left co-ordinate of the cell
			 */
	size_t row,      /* Row number (relates to y position). */
	size_t col,      /* Column number (relates to x position). */
	double *x,      /* write-only. Returns x of true co-ordinate */
	double *y)      /* write-only. Returns y of true co-ordinate */
{
	return RrowCol2Coords(m,
		(double)row+(inCellPos ? 0.5 : 0.0),
		(double)col+(inCellPos ? 0.5 : 0.0),
		x,y);
}

/* compute true world co-ordinate from row, column index
 * RrowCol2Coords computes the true world co-ordinate from a 
 * row, column index. The row,column co-ordinate can be fractions.
 * For example (row,col) = (0.5,0.5) computes the (x,y) co-ordinate of
 * the centre of the upper left pixel. Secondly, the row and column co-ordinate
 * don't have to be on the map. They are just relative to upper left position.
 * For example (row,col) = (-0.5,0.5) computes the (x,y) co-ordinate of
 * the centre of the pixel that is right above upper left pixel. 
 *
 * returns
 *  0  if the co-ordinate is outside the map.
 *  1 if inside.
 * -1 in case of an error.
 *
 * Merrno
 * ILL_CELLSIZE
 */
int RrowCol2Coords(const MAP *m, /* map handle */
		    double row,  /* Row number (relates to y position). */
		    double col,  /* Column number (relates to x position). */
		    double *x,   /* write-only. x co-ordinate */
		    double *y)   /* write-only. y co-ordinate */
{
	double cs     = m->raster.cellSizeX;
	double c      = m->angleCos;
	double s      = m->angleSin;
	double yRow   = cs * row;
	double xCol   = cs * col;
	double xCol_t = xCol * c - yRow * s;
	double yRow_t = xCol * s + yRow * c;

	if (m->raster.cellSizeX <= 0 || m->raster.cellSizeY <= 0
	    || m->raster.cellSizeX != m->raster.cellSizeX )
	{
		M_ERROR(ILL_CELLSIZE);
		goto error;
	}

	*x = m->raster.xUL + xCol_t;
	if (m->main.projection == PT_YINCT2B)
		*y = m->raster.yUL + yRow_t;
	else  /* all other projections */
		*y = m->raster.yUL - yRow_t;

	return( (m->raster.nrRows > row) && (m->raster.nrCols > col) &&
	        (row >= 0) && (col >= 0));
error:  return(-1);
}
