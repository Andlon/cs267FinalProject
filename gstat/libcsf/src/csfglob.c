/*
 * csfglob.c
 */


#include "csf.h"
#include "csfimpl.h"

#ifndef lint
 static const char *rcs_id = 
 "$Header: /home/pcrtree/SRC.RCS/libs/csf/RCS/csfglob.c,v 2.1 1996/12/29 19:35:21 cees Exp $";
 static const char *includeRcsId[4] = { RCS_ID_CSF_H, RCS_ID_CSFIMPL_H, 
                                 RCS_ID_CSFLEGEN_H , RCS_ID_CSFTYPES_H };
#endif

/* global variable set on the last error condition
 * Most functions sets this variable in case of an error condition.
 */
int Merrno = NOERROR; 
