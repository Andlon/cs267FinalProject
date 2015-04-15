/*
 * rextend.c
 */
#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/rextend.c,v 1.4 1997/10/07 09:45:40 cees Exp $";
#endif

#include "csf.h"
#include "csfimpl.h"
#include <math.h>

static double FmodRound(
	double v,
	double round)
{
	double rVal = fmod(v, round);
	if(rVal == 0)
		return round;
	else
		return rVal;
}

/* compute (xUL,yUL) and nrRows, nrCols from some coordinates
 * RcomputeExtend computes parameters to create a raster maps
 * from minimum and maximum x and y coordinates, projection information,
 * cellsize and units. The resulting parameters are computed that the 
 * smallest raster map can be created that will include the two 
 * coordinates given, assuming a default angle of 0.
 * Which coordinates are the maximum or minimum are
 * determined by the function itself.
 */
void RcomputeExtend(
	REAL8 *xUL,     /* write-only, resulting xUL */
	REAL8 *yUL,     /* write-only, resulting yUL */
	size_t *nrRows, /* write-only, resulting nrRows */
	size_t *nrCols, /* write-only, resulting nrCols */
	double x_1,      /* first x-coordinate */ 
	double y_1,      /* first y-coordinate */
	double x_2,      /* second x-coordinate */
	double y_2,      /* second y-coordinate */
	CSF_PT projection, /* required projection */
	REAL8 cellSize, /* required cellsize, > 0 */
	double rounding) /* assure that (xUL/rounding), (yUL/rouding)
	                  * (xLL/rounding) and (yLL/rounding) will
	                 * will all be an integers values > 0 
	                 */
{
    /*
     * xUL ______
	   |    |
	   |    |
	   |    |
	   ------

     */
	double yLL,xUR = x_1 > x_2 ? x_1 : x_2;
	*xUL = x_1 < x_2 ? x_1 : x_2;
	*xUL -= FmodRound(*xUL, rounding);
	xUR  += rounding - fmod(xUR, rounding);
	POSTCOND(*xUL <= xUR);
	*nrCols = (size_t)ceil((xUR - *xUL)/cellSize);
	if (projection == PT_YINCT2B)
	{
		 yLL = y_1 > y_2 ? y_1 : y_2;  /* highest value at bottom */
		*yUL = y_1 < y_2 ? y_1 : y_2;  /* lowest value at top */
	        *yUL -= FmodRound(*yUL, rounding);
	         yLL += rounding-fmod( yLL, rounding);
	}
	else
	{
		 yLL = y_1 < y_2 ? y_1 : y_2;  /* lowest value at bottom */
		*yUL = y_1 > y_2 ? y_1 : y_2;  /* highest value at top */
	        *yUL += rounding-fmod(*yUL, rounding);
	         yLL -= FmodRound( yLL, rounding);
	}
	*nrRows = (size_t)ceil(fabs(yLL - *yUL)/cellSize);
}
