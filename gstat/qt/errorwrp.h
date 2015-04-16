#if defined(__cplusplus)
extern "C" {
#endif

#include <setjmp.h>

void gstat_error_jmp(const char *mess, int level);
extern jmp_buf gstat_buf;

#if defined(__cplusplus)
}
#endif
