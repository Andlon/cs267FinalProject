#ifdef __cplusplus 
extern "C" {
#endif

#ifndef INCLUDED_APP_GSTAT

/* dummy app_GStat struct decl: */
typedef struct {
	int xxx; /* BORLANDC gives error on empty struct */
} app_GStat;

#endif

extern app_GStat *pt;

extern void setProgress(app_GStat *p, int n, int total);
extern void showWarning(app_GStat *p, const char *str);
extern void logMessage(app_GStat *p, const char *str);
extern void set_pointer(app_GStat *p);

extern void wrap_progress(int n, int total);
extern void wrap_warning(const char *str);
extern void wrap_log(const char *str);

#ifdef __cplusplus 
}
#endif
