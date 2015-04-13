/* src/config.h.  Generated from config.hin by configure.  */
/* config.hin.  Generated automatically from configure.in by autoheader.  */
#ifndef CONFIG_H
# define CONFIG_H

/* Define to empty if the keyword does not work.  */
/* #undef const */

/* Define USING_R if compiled as an R package */
/* #undef USING_R */

/* Define if you don't have vprintf but do have _doprnt.  */
/* #undef HAVE_DOPRNT */

/* Define if you have the vprintf function.  */
#define HAVE_VPRINTF 1

/* Define if you have the snprintf function.  */
#define HAVE_SNPRINTF 1

/* Define if you have the vsnprintf function.  */
#define HAVE_VSNPRINTF 1

/* Define if ansi-c headers are present. */
#define STDC_HEADERS 1

/* Define if both time.h and sys/time.h can be included */
#define TIME_WITH_SYS_TIME 1

/* Define if sys/time.h is present */
/* #undef HAVE_SYS_TIME */

/* Define if you have the HAVE_<xxx> function.  */
#define HAVE_STRSTR 1
#define HAVE_STRTOD 1
#define HAVE_STRTOL 1
#define HAVE_POPEN 1
#define HAVE_GETTIMEOFDAY 1
#define HAVE_DRAND48 1
#define HAVE_GETOPT 1
#define HAVE_GETOPT_LONG 1

/* Define if you have the <unistd.h> header file.  */
#define HAVE_UNISTD_H 1

/* Define if you have the <unistd.h> header file.  */
#define HAVE_GETOPT_H 1

/* Define if you have the GNU scientific library (-lgsl).  */
/* #undef HAVE_LIBGSL */

/* for the PCRASTER cross system format (csf) library (-lcsf).  */
/* #undef HAVE_LIBCSF */

/* Define if you have the ncurses library (-lncurses).  */
#define HAVE_LIBNCURSES 1

/* Define if you have the curses library (-lcurses).  */
/* #undef HAVE_LIBCURSES */

/* Define if you have the gnu readline library (-lreadline).  */
/* #undef HAVE_LIBREADLINE */

/* Define if you have the gd library (-lgd); */
/* #undef HAVE_LIBGD */

/* Define if you have the gdal library (-lgdal); */
#define HAVE_LIBGDAL 1

/* ... with PNG support...  */
/* #undef HAVE_GDIMAGEPNG */

/* ...or with GIF support. */
/* #undef HAVE_GDIMAGEGIF */

/* Define if you have the netcdf library (-lnetcdf) for GMT map support  */
#define HAVE_LIBNETCDF 1

/* Define if you want to link/include extdbase extension */
/* #undef HAVE_EXT_DBASE */

/* Define if words are big endian */
/* #undef WORDS_BIGENDIAN */

/* include sparse matrix routines in meschach */
/* #undef HAVE_SPARSE */

/* define if target is libgstat.a */
/* #undef LIBGSTAT */

/* define if GRASS' libgis.a is present */
/* #undef HAVE_LIBGIS */

#endif
