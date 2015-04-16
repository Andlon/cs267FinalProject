#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cgic.h"

#define GSTAT "/www/gstat/bin/gstat"
#define GNUPLOT "/www/gstat/bin/gnuplot"

#define STREQ(a,b) (strcmp(a,b) == 0)

int cgiMain() {
	double b, B, s, S, nug, sill, range;
	char vgm[100], *models[] = { "Exp", "Sph", "Gau", "Cir" },
		call[100], args[100], output_value[10];
	int whichmodel = 0, exitstatus, nbh, out = 0;

#if DEBUG
	/* Load a saved CGI scenario if we're debugging */
	cgiReadEnvironment("/home/boutell/public_html/capcgi.dat");
#endif
	/* parse input: */
	cgiFormDoubleBounded("mins", &s, 1.0e-10, 1.0e38, 1);
	cgiFormDoubleBounded("maxs", &S, 1.0e-10, 1.0e38, 10);
	cgiFormDoubleBounded("minb", &b, 0.0, 1.0e38, 1);
	cgiFormDoubleBounded("maxb", &B, 1.0e-10, 1.0e38, 10);
	cgiFormDoubleBounded("nug", &nug, 0.0, 1.0e38, 1);
	cgiFormDoubleBounded("sill", &sill, 0.0, 1.0e38, 1);
	cgiFormDoubleBounded("range", &range, 1e-10, 1.0e38, 10);
	cgiFormIntegerBounded("nbh", &nbh, 1, 100, 25);
	cgiFormSelectSingle("model", models, 4, &whichmodel, 0);
	cgiFormString("out", output_value, 10);
	sprintf(vgm, "%g Nug() + %g %s(%g)", nug, sill, models[whichmodel], range);
	if (strcmp(output_value, "graph") == 0) {
		cgiHeaderContentType("image/gif");
		fflush(cgiOut);
		sprintf(call, 
				"%s -e ossfim -m '-' -s %g -S %g -b %g -B %g -n %d -v '%s'\n", 
				GSTAT, s, S, b, B, nbh, vgm);
		system(call); /* outputs a GIF only */
	} else if (strcmp(output_value, "vgm") == 0) {
		cgiHeaderContentType("image/gif");
		fflush(cgiOut);
		sprintf(call, "%s -e ossfim -V '%s' | %s\n", GSTAT, vgm, GNUPLOT);
		system(call);
	} else {
		cgiHeaderContentType("text/html");
		fprintf(cgiOut, "<HTML><HEAD>\n");
		fprintf(cgiOut, "<TITLE>gstat output</TITLE></HEAD>\n");
		fprintf(cgiOut, "<BODY>\n");

		/* plot graph? */
		if (STREQ(output_value, "all") || STREQ(output_value, "outgraph")) {
			/* graph: */
			sprintf(call, "ossfimcgi?mins=%g&maxs=%g&minb=%g&maxb=%g&nug=%g&sill=%g&model=%s&range=%g&nbh=%d&out=graph",
			s, S, b, B, nug, sill, models[whichmodel], range, nbh);
			fprintf(cgiOut, "<p><h2>Graph</h2> (x = sample spacing; y = block size; colour denotes block kriging standard (prediction) error):\n");
			fprintf(cgiOut, "<p><img src=\"%s\">\n", call);
		}
		if (STREQ(output_value, "all") || STREQ(output_value, "outvgm")) {
			/* variogram: */
			sprintf(call, "ossfimcgi?mins=%g&maxs=%g&minb=%g&maxb=%g&nug=%g&sill=%g&model=%s&range=%g&nbh=%d&out=vgm",
			s, S, b, B, nug, sill, models[whichmodel], range, nbh);
			fprintf(cgiOut, "<p><h2>Variogram:</h2>\n");
			fprintf(cgiOut, "<p><img src=\"%s\">\n", call);
		}
		if (STREQ(output_value, "all") || STREQ(output_value, "outdata")) {
			fprintf(cgiOut, "<h2>gstat output:</h2>\n");
			/* call gstat: */
			sprintf(args, "-e ossfim -s %g -S %g -b %g -B %g -n %d -v '%s'\n", 
				s, S, b, B, nbh, vgm);
			fprintf(cgiOut, "gstat call:<p><pre>gstat %s\n</pre>", args);
			fprintf(cgiOut, "<p>[column 1 = sample spacing; column 2 = block size; column 3 = block kriging prediction error]:<p><pre>\n");
			fflush(cgiOut);
			sprintf(call, "%s %s", GSTAT, args);
			exitstatus = system(call);
			fprintf(cgiOut, "</pre>\n");
			if (exitstatus)
				fprintf(cgiOut, "exit status: %d\n", exitstatus);
		}
		fprintf(cgiOut, "<p>ossfimcgi uses <a href=\"http://www.boutell.com/cgic/\">cgic</a>.\n");
		fprintf(cgiOut, "</BODY></HTML>\n");
	}
	return 0;
}
