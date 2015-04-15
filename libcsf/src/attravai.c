/*
 * attravai.c
 */
#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/attravai.c,v 2.1 1996/12/29 19:35:21 cees Exp $";
#endif

#include "csf.h"
#include "csfimpl.h"

/* check if an attribute is available
 * MattributeAvail search for the given id in the map.
 *
 * returns
 *  0 if the attribute is not available,
 *  nonzero if the attribute is available
 *
 * Merrno
 * ILLHANDLE
 */
int MattributeAvail(
	MAP *m,    /* map handle */
	CSF_ATTR_ID id)  /* identification of attribute */
{
	ATTR_CNTRL_BLOCK b;

     	if (! CsfIsValidMap(m))
	 	return 0;
	return(CsfGetAttrBlock(m, id, &b) != 0);
}
