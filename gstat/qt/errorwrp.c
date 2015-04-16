#include "errorwrp.h"

jmp_buf gstat_buf;

void gstat_error_jmp(const char *mess, int level) {
	/* jump back to where the jump buffer was set: */
	longjmp(gstat_buf, level == 0 ? -1 : level);
}
