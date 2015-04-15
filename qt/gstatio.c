#include <assert.h>
#include <stdio.h>
#include "gstatio.h"

static app_GStat *pt = NULL;

void set_pointer(app_GStat *p) {
	assert(p);
	pt = p;
}

void wrap_progress(int n, int total) {
	assert(pt);

	setProgress(pt, n, total);
}

void wrap_warning(const char *str) {
	assert(pt);

	showWarning(pt, str);
}

void wrap_log(const char *str) {
	assert(pt);

	logMessage(pt, str);
}
