/************************************************************************/
/*                                                                      */
/* Things that should be in every module                                */
/*   Only one time include                                              */
/*                                                                      */
/*                                                                      */
/************************************************************************/
#ifndef __STDDEFX
#define __STDDEFX

#ifdef INTEL32
#include <unistd.h>
#endif /* INTEL32 */

/* ANSI */

#ifndef ANSI
#define ANSI
#endif

#ifdef ANSI
#include <stddef.h>
#include <stdlib.h>
/* #include <float.h> */
#include <limits.h>

#include <errno.h>
#endif /* ANSI */

/* SELF */

#include "environ.h"
#include "typedef.h"
#include "debug.h"
/*#include "stdlibx.h" */	/*extentions*/

/***********************************************/
/* Extended library functions                  */
/***********************************************/

#ifndef MIN
#define MIN(a,b)	((a < b) ? a : b)
#endif

#ifndef MAX
#define MAX(a,b)	((a > b) ? a : b)
#endif

#ifndef MAX
#define ABS(a) 		((a >= 0) ? (a) : (-a))
#endif

#define ERRNO(x)	(errno == x)

#define TABLE_INDEX(table, type, i)	((type *)*(VOID **)table)[i]

#ifdef DIAGNOSE
#define IFDIAGNOSE(action)     	action
#else
#define IFDIAGNOSE(action)
#endif

#endif /*__STDDEFX*/
