/*
 * rrowcol.c
 */
#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/rrowcol.c,v 2.1 1996/12/29 19:35:21 cees Exp $";
#endif

#include "csf.h"
#include "csfimpl.h"

#include <math.h> /* floor */

/* compute (fractional) row, column index from true world co-ordinate.
 * Rcoord2RowCol computes row, column index from true world co-ordinate. 
 * The row and column co-ordinate are returned as fractions (See parameters
 * section).
 * The x,y co-ordinate
 * don't have to be on the map. They are just relative to upper left position.
 *
 * returns
 * 0 if the co-ordinate is outside the map,
 * 1 if inside,
 * -1 in case of an error
 *
 * Merrno
 * ILL_CELLSIZE
 */
int Rcoords2RowCol(
 const MAP *m,  /* map handle */
 double x,      /* x of true co-ordinate */
 double y,      /* y of true co-ordinate */
 double *row,   /* write-only. Row index (y-pos). floor(row) is row number,
                 * if (row >= 0) then fmod(row, 1) is in-pixel displacement from pixel-top,
                 * if (row <0) then fmod(row, 1) is in-pixel displacement from pixel-bottom. 
                 */
 double *col)   /* write-only. Column index (x-pos). floor(col) is column number,
                 * if (col >= 0) then fmod(col, 1) is in-pixel displacement from pixel-left, 
                 * if (col <0) then fmod(col, 1) is in-pixel displacement from pixel-right.
                 */
{
	double row_,col_; /* use copies, func( , , ,&dummy, &dummy) will fail
	                   * otherwise
	                   */
	double cs = m->raster.cellSizeX;
	double xCol = (x - m->raster.xUL) / cs; 
	double yRow = (((m->main.projection == PT_YINCT2B) 
	                ? (y - m->raster.yUL)
	                : (m->raster.yUL - y)) / cs); 
	/* rotate clockwise: */
	double c = m->angleCos;    /* cos(t) == cos(-t) */
	double s = -(m->angleSin); /* -sin(t) == sin(-t) */
	col_  = xCol * c - yRow * s;
	row_  = xCol * s + yRow * c;

	if (m->raster.cellSizeX <= 0 || m->raster.cellSizeY <= 0
	    || (m->raster.cellSizeX != m->raster.cellSizeY ) )
	{ /* CW we should put this in Mopen */ 
		M_ERROR(ILL_CELLSIZE);
		goto error;
	}
        *row = row_;
        *col = col_;
	return( row_ >= 0 && col_ >= 0 && 
	        (m->raster.nrRows > row_) && (m->raster.nrCols > col_) );
error:  return(-1);
}


/* compute row, column number of true world co-ordinate
 * RgetRowCol computes row, column number of true world co-ordinate.
 *
 * returns
 * 0  if the co-ordinate is outside the map,
 * 1 if inside,
 * -1 in case of an error
 *
 * Merrno
 * ILL_CELLSIZE
 */
int RgetRowCol(
	const MAP *m, /* map handle */
	double x,     /* x of true co-ordinate */
	double y,     /* y of true co-ordinate */
	size_t *row,   /* write-only. Row number (y-pos).
	               * Undefined if (x,y) is outside of map
	               */
	size_t *col)   /* write-only. Column number (x-pos).
	               * Undefined if (x,y) is outside of map
	               */
{
	double row_d,col_d;
	int    result;
	result = Rcoords2RowCol(m,x,y,&row_d,&col_d);
	if (result > 0)
	{
		*row = (size_t)floor(row_d);
		*col = (size_t)floor(col_d);
	}
	return(result);
}
