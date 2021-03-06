/*
 * attrsize.c
 */
#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/attrsize.c,v 2.1 1996/12/29 19:35:21 cees Exp $";
#endif

#include "csf.h"
#include "csfimpl.h"

/* get the size of an attribute (LIBRARY_INTERNAL)
 * returns
 * 0 if the attribute is not available,
 * or the nonzero size if the attribute is available.
 */
size_t CsfAttributeSize(
	 MAP   *m,    /* map handle */
	 CSF_ATTR_ID id)    /* identification of attribute */
{
	ATTR_CNTRL_BLOCK b;

	if (CsfGetAttrBlock(m, id, &b) != 0)
		return b.attrs[CsfGetAttrIndex(id, &b)].attrSize;
        return 0;
}
