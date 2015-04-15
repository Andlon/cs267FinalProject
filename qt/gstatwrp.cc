/* 
 * gui.cc: interface from gstat to a GUI. Started in May 1998.
 */

/*! \file gstatwrp.cc
 *  \brief gstat C++ functions wrappers
 *
 * This module provides an abstract gui interface layer to gstat:
 * it gives access to all that's needed from gstat for a gui, 
 * but it doesn't do any graphics or events.
 *
 * We will implement a gstat class, see the corresponding .h file
 * for a definition.
 */

/*
 * TODO: setup_minmax() and calc_polynomials()
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <setjmp.h>

#include "defs.h"
#include "debug.h"
#include "lex.h"
#include "utils.h"
#include "userio.h"
#include "data.h"
#include "vario.h"
#include "sem.h"
#include "fit.h"
#include "glvars.h"
#include "version.h"
#include "writecmd.h"
#include "errorwrp.h"

#include "gstatwrp.h"
#include "gstatexcept.h" // throw

#define setup_error_handler(); \
	{ \
	int error_level; \
	reset_gstat_errno(); \
	if ((error_level = setjmp(gstat_buf)) == 0) \
		set_gstat_error_handler(gstat_error_jmp);  \
	else  \
		process_gstat_error(error_level); \
	}

int Gstat::active = 0;

Gstat::Gstat(app_GStat *p) 
{
 	assert(!active);
	assert(p);

 	active = 1;

	init_userio(0);
	init_global_variables();
	// logfile_name = "gstat.log";
	// parse_gstatrc();
	gl_register_pairs = 0; // would dependent on data list length...
	set_method(UIF); // suppress printing in ?print_cmd()

	// set call-back functions for pr_warning, printlog and print_progress
	set_pointer(p);
	set_gstat_warning_handler(wrap_warning);
	set_gstat_log_handler(wrap_log);
	set_gstat_progress_handler(wrap_progress);
}

Gstat::~Gstat() 
{
	remove_all(); // clean up the ``private'' data in glvars.c
	active = 0; // enable restart
	return;
}

const char *Gstat::get_version() const
{
	return VERSION;
}

void Gstat::process_gstat_error(int level)
{
	switch (level) {
		case ER_NULL:
		case ER_IMPOSVAL:
		case ER_RANGE:
		case ER_MESCHACH:
		case ER_VARNOTSET:
			throw GStatFatalError(get_gstat_error_message());
			break;
		case ER_SYNTAX:
		case ER_NOFILE:
		case ER_WRITE:
		case ER_READ:
		case ER_RDFLT:
		case ER_RDINT:
		case ER_ARGOPT:
		case ER_DOMAIN:
		case ER_MEMORY:
		case ER_SECURE:
		case ER_PWRITE:
		case ER_PREAD:
			throw GStatNormalError(get_gstat_error_message());
			break;
		case ER_NOERROR:
		case ER_NOCMD:
		case ER_IO:
		case ER_NOCURSES:
		default:
			assert(0);
	}
}

/*!
 * add a data variable
 * \return always returns 0 (or breaks on an error)
 */
void Gstat::add_data(const char *name, const char *fname, int x, int y, int v) 
{ 
	DATA **d;
	int id;

	assert(name);
	assert(fname);

	id = which_identifier(name); // resizes data if necessary
	d = get_data();
	d[id]->id = id; /* the magic phrase... */

	d[id]->fname = string_dup(fname);
	d[id]->colnx = x;
	d[id]->colny = y;
	d[id]->colnvalue = v;

	setup_error_handler();

	read_data(d[id]);

	d_set(name);

	return;
}

void Gstat::add_data(int id) 
{
	DATA **d;

	assert(id >= 0);
	assert(id < get_n_vars()); // NOT the member function get_nvars()

	d = get_data();
	d[id]->id = id; /* the magic phrase... */

	read_data(d[id]);
	d_set(id);
	return;
}

/*
 * remove variable name 
 */
void Gstat::remove_data(const char *name) 
{
	int id;

	assert(name);
	assert(get_nvars() > 0);
	id = which_identifier(name);
	if (dptr->id == id) {
		if (get_nvars() > 1) // set dptr one back
			d_set(id - 1); // don't dangle!
		else
			dptr = NULL; // this may be helpful later on.
	}
	remove_id(id); // releases memory of id
	return;
}

/* 
 * report summary statistics of data variable to file f
 */
void Gstat::report(FILE *f, const char *name) const
{
	assert(name);
	report(f, which_identifier(name));
}

void Gstat::report(FILE *f, int id) const
{
	DATA **d;

	assert(id >= 0);
	assert(id < get_n_vars());

	d = get_data();
	fprintf(f, "data(%s): ", name_identifier(id));
	report_data(d[id]);
}

size_t Gstat::get_nvars() const
{
	return get_n_vars(); /* gstat C function */
}

//! \return 1 on error, 0 on success
int Gstat::read_commands(const char *cmd, const char *filename) 
{
	setup_error_handler();

	assert(get_n_vars() == 0);

	if (cmd == NULL) {
		assert(filename != NULL);
		cmd = string_file(filename);
	} else {
		if (filename == NULL)
			filename = "?";
	}

	if (parse_cmd(cmd, filename))
		ErrMsg(ER_SYNTAX, filename);	

	if (get_n_vars() > 1) { // remove all but the first:
		pr_warning("multiple data definitions are ignored");
		while(get_n_vars() > 1)
			remove_id(get_n_vars() - 1);
	}

	if (get_n_vars() > 0) {
		VARIOGRAM *v = get_vgm(LTI(0,0));
		if (v->n_models > 0) {
			v->n_models = 0;
			pr_warning("Variogram definition(s) are ignored");
		}
	} else
		pr_warning("%s contained no data definition", filename);

	// assert(get_n_vars() == 1);

	// in case one or more data(xx) commands appeared in cmd:
	// (will effectively only read one data var at this moment)
	for (int i = 0; i < get_n_vars(); i++)
		add_data(i);

	return 0;
}

void Gstat::write_commands(const char *fname) //!< write command file
{
	FILE *f;

	assert(fname != NULL);
	setup_error_handler();

	f = efopen(fname, "w");
	fprint_cmd(f);
	efclose(f);
}

void Gstat::d_set(const char *name)
{
	int id;

	assert(name);
	id = which_identifier(name);
	d_set(id);
	return;
}

void Gstat::d_set(int id)
{
	DATA **d;
	d = get_data();
	dptr = d[id];
}

const char *Gstat::d_name(int id) const
{
	assert(id >= 0);
	assert(id < get_n_vars());
	return name_identifier(id);
}

size_t Gstat::d_id(const char *name) const
{
	assert(name);
	return (size_t) which_identifier(name);
}

size_t Gstat::d_get_n() const // number of data points in this variable
{
	assert(dptr);
	return (size_t) dptr->n_list;
}

void Gstat::d_get_xyv(int i, double *x, double *y, double *v) const
{
	assert(dptr);
	assert(dptr->id >= 0);
	assert(i >= 0);
	assert(i < dptr->n_list);

	*x = dptr->list[i]->x;
	*y = dptr->list[i]->y;
	*v = dptr->list[i]->attr;
	return;
}

void Gstat::v_set(const char *name) {
	v_set(name, name);
}

void Gstat::v_set(const char *a, const char *b) {
	int id1, id2;

	if (a == b)
		id1 = id2 = which_identifier(a);
	else {
		id1 = which_identifier(a);
		id2 = which_identifier(b);
	}
	vptr = get_vgm(LTI(id1,id2));
	vptr->ev->evt = (id1 == id2 ? SEMIVARIOGRAM : CROSSVARIOGRAM);
	vptr->id = LTI(id1, id2);
	vptr->id1 = id1;
	vptr->id2 = id2;
}

void Gstat::v_set_cutoff(double val)
{
	assert(val > 0.0);
	assert(vptr->ev);

	vptr->ev->cutoff = gl_cutoff = val;
	vptr->ev->recalc = vptr->ev->refit = 1;
}

double Gstat::v_get_cutoff() const
{
	assert(vptr->ev);
	return vptr->ev->cutoff;
}

void Gstat::v_set_width(double val)
{
	assert(val >= 0.0);
	assert(vptr->ev);

	vptr->ev->iwidth = gl_iwidth = val;
	vptr->ev->recalc = vptr->ev->refit = 1;
}

double Gstat::v_get_width() const
{
	assert(vptr->ev);
	return vptr->ev->iwidth;
}

void Gstat::v_calc(const char *name) 
{ // calculate variogram of variable name
	v_calc(name, name);
}

void Gstat::v_calc(const char *a, const char *b)
{ // calculate cross variogram of variables a and b 

	assert(a);
	assert(b);
	v_set(a,b);
	assert(vptr);

	setup_error_handler();

	calc_variogram(vptr, NULL); /* use "" to write to stdout */
}

void Gstat::v_fit()
{
	FIT_TYPE ft;

	assert(vm_has_model());

	setup_error_handler();

	setup_meschach_error_handler();

	ft = (FIT_TYPE) gl_fit;
	if (ft > NO_FIT && ft <= MIVQUE_FIT)
		vptr->ev->fit = ft;
	else
		vptr->ev->fit = WLS_FIT_MOD;

	vptr->ev->refit = 1;
	fit_variogram(vptr);
	update_variogram(vptr);
}

double Gstat::v_get_model(double dx, double dy, double dz) const
{ // get model value for vector (dx,dy,dz)
	assert(vptr);
	assert(vptr->n_models > 0);
	return get_semivariance(vptr, dx, dy, dz);
}

size_t Gstat::v_get_n_est() const
{ // get number of sample variogram estimates
	assert(vptr);
	assert(vptr->ev);
	return((size_t) vptr->ev->n_est);
}

double Gstat::v_get_est(size_t i) const
{ // get sample variogram value i
	assert(vptr);
	assert(vptr->ev);
	assert(i < (size_t) vptr->ev->n_est);
	return vptr->ev->gamma[i];
}

double Gstat::v_get_dist(size_t i) const
{ // get average distance for value i
	assert(vptr);
	assert(vptr->ev);
	assert(i < (size_t) vptr->ev->n_est);
	return vptr->ev->dist[i];
}

size_t Gstat::v_get_npp(size_t i) const
{ // get number of pairs of points for value i
	assert(vptr);
	assert(vptr->ev);
	assert(i < (size_t) vptr->ev->n_est);
	return vptr->ev->cloud ? 1 : (size_t) vptr->ev->nh[i];
}

void Gstat::v_get_ppair(size_t i, size_t j,
		double *x1, double *y1, double *x2, double *y2) const
{ // for j = 0 ... get_n_pairs_i(i), point pairs (x1,y1)--(x2,y2) 
	DPOINT *dpt1, *dpt2;
	DATA **d, *d1, *d2;

	assert(vptr);
	assert(vptr->ev);
	assert(i < (size_t) vptr->ev->n_est);
	assert(j < vptr->ev->nh[i]);
	assert(gl_register_pairs);
	// assert(!vptr->ev->cloud);
	// assert(vptr->id1 == vptr->id2); 
	if (! vptr->ev->cloud)
	{
		dpt1 = vptr->ev->pairs[i][j*2];
		dpt2 = vptr->ev->pairs[i][(j*2)+1];
	}
	else
	{
		assert(j == 0);
		d = get_data();
		d1 = d[vptr->id1];
		d2 = d[vptr->id2];
		dpt1 = d1->list[LOW_NH(vptr->ev->nh[i])];
		dpt2 = d2->list[HIGH_NH(vptr->ev->nh[i])];
	}
	*x1 = dpt1->x;
	*y1 = dpt1->y;
	*x2 = dpt2->x;
	*y2 = dpt2->y;
}

// vm_* function to query/set variogram model stuff
int    Gstat::vm_has_model() const // non-zero if a model is present
{
	return(vptr->n_models > 0);
}

const char **Gstat::vm_get_possible_models() const // NULL-terminated char* array
{
	int i, j, n;
	const char **cpp;

	for (i = n = 1; v_models[i].model != NOT_SP; i++)
		n++;

	cpp = (const char **) emalloc(n * sizeof(const char *));

	for (i = 1, j = 0; i < n; i++) {
		switch (v_models[i].model) { /* provide a ``save'' subset */
			case EXPONENTIAL: 
			case GAUSSIAN: 
			case BESSEL:
			case PENTASPHERICAL:
			case SPHERICAL:
			case CIRCULAR:
			case POWER:
			/* case LINEAR: -- solve the Lin(0) case first */
				cpp[j] = v_models[i].name_long;
				j++;
				break;
			default: 
				break;
		}
	}
	cpp[j] = NULL;
	return cpp;
}

const char *Gstat::vm_get_nugget_string() const
{
	int i;

	for (i = 0; v_models[i].name_long != NULL; i++)
		if (v_models[i].model == NUGGET)
			return v_models[i].name_long;
	assert(0);
	return "";
}

const char *Gstat::vm_get_model_description() const
{
	if (vptr->descr != NULL)
		return vptr->descr;
	else
		return "";
}

int    Gstat::vm_add_model(const char *name, double sill, double range) 
{
	return add_edit_model(-1, name, sill, range);
}

int    Gstat::add_edit_model(int id, const char *name, double sill, double range)
{

	VGM_MODEL p;
	int retval;

	init_variogram_part(&p);
	p.model = which_variogram_model(name);

	assert(p.model != NOT_SP);

	p.sill = sill;
	p.range = range;
	p.id = id;
	retval = push_variogram_model(vptr, p); // will search for p.id in list
	update_variogram(vptr);
	return retval;
}

void   Gstat::vm_edit_model(int id, const char *name, double sill, double range)
{
	add_edit_model(id, name, sill, range);
	return;
}

void   Gstat::vm_delete_model(int id) 
{
	int i, where = -1;

	for (i = 0; i < vptr->n_models; i++) {
		if (vptr->part[i].id == id) {
			where = i;
			break;
		}
	}

	assert(where >= 0);

	for (i = where; i < vptr->n_models - 1; i++)
		vptr->part[i] = vptr->part[i+1];
	vptr->n_models--;
	update_variogram(vptr);
}

double Gstat::vm_get_default_range(const char *model) const
{
	size_t i, n;
	double mean = 0.0, factor;

	VGM_MODEL_TYPE m = which_variogram_model(model);

	assert(m != NOT_SP);

	switch (m) {
		case EXPONENTIAL: 
			factor = 3.0; 
			break;
		case GAUSSIAN: 
			factor = sqrt(3.0); 
			break;
		case BESSEL: 
			factor = 4.0; 
			break;
		case POWER: 
			return 1.0; /* EXIT() function */
		default: 
			factor = 1.0; 
			break;
	}

	if (v_get_n_est() > 0) {
		for (i = n = 0; i < v_get_n_est(); i++) {
			if (v_get_npp(i)) {
				mean += v_get_dist(i);
				// printf("(%g)\n", v_get_dist(i));
				n++;
			}
		}
		return mean / (factor * n);
	}
	return 1.0;
}

double Gstat::vm_get_sill(int id) const
{
	int i;
	for (i = 0; i < vptr->n_models; i++)
		if (id == vptr->part[i].id)
			return vptr->part[i].sill;
	assert(0);
}

double Gstat::vm_get_range(int id) const
{
	int i;
	for (i = 0; i < vptr->n_models; i++)
		if (id == vptr->part[i].id)
			return vptr->part[i].range;
	assert(0);
}

double Gstat::vm_get_default_sill() const
{
	size_t i, n;
	double mean = 0.0;

	if (v_get_n_est() > 0) {
		for (i = n = 0; i < v_get_n_est(); i++) {
			if (v_get_npp(i)) {
				mean += v_get_est(i);
				n++;
			}
		}
		return 0.75 * mean / n;
	}
	return 1.0;
}

double Gstat::vm_get_default_nugget() const 
{
	return 0.5 * vm_get_default_sill();
}

void Gstat::vm_print() const
{
	printf("[%s]\n", vptr->descr);
}

void Gstat::set_warning_handler(void (*warning_fn)(const char *message)) 
{
	set_gstat_warning_handler(warning_fn);
}

void Gstat::set_error_handler(void (*error_fn)(const char *message, int level)) 
{
	set_gstat_error_handler(error_fn);
}

#ifdef MAIN_WRP
/*
 * main program to test Gstat member functions
 * make with `make gui'
 */
int main(int argc, char *argv[]) 
{
	size_t i, j; //, id1, id2, id3;

	for (i = 1; i < (size_t) argc; i++) {
		try {
			Gstat g(argc, argv);
			printf("(%s:)\n", argv[i]);
			if (g.read_commands(string_file(argv[i]), argv[i]))
				printf("an error occured!\n");
			for (j = 0; j < g.get_nvars(); j++) 
			{
				printf("[[%s]]:\n", g.d_name(j));
				g.report(stdout, j);
				g.v_set(g.d_name(j));
				g.v_calc(g.d_name(j));
				printf("[%d] ", g.v_get_n_est());
				// for (size_t k = 0; k < g.v_get_n_est(); k++) 
				// {
					// for (size_t l = 0; l < g.v_get_npp(k); l++) 
					// {
					// double x1, x2, y1, y2;
					// g.v_get_ppair(k,l,&x1,&y1,&x2,&y2);
					// printf("[%d-%d:%g,%g-%g,%g]", k,l,x1,y1,x2,y2);
					// }
				// }
				printf("\n");
			}
			printf("[ok ... ]\n");
			while (g.get_nvars())
				g.remove_data(g.d_name(0));
	
		}

		catch (GStatFatalError) {
			printf("caught a GstatFatalError!\n");
		}
		catch (GStatNormalError) {
			printf("caught a GstatNormalError!\n");
		}
		catch (...) {
			printf("caught an (...) error!\n");
		}
	}
	
	try {
		Gstat g(argc, argv);
		g.add_data("xyz", "zinc.eas", 1, 2, 3);
		g.v_set("xyz");
		g.v_calc("xyz");
		// g.v_set_width(1);
		// g.v_set_cutoff(100);
		g.v_calc("xyz");
		printf("cutoff: %g\n", g.v_get_cutoff());
		for (i = 0; i < g.v_get_n_est(); i++)
			if (g.v_get_npp(i))
				printf("%g %g\n", g.v_get_dist(i), g.v_get_est(i));
		g.v_set_width(20);
		g.v_set_cutoff(200);
		g.v_calc("xyz");
		printf("cutoff: %g\n", g.v_get_cutoff());
		for (i = 0; i < g.v_get_n_est(); i++)
			printf("%g %g %d\n", g.v_get_dist(i), g.v_get_est(i),
				g.v_get_npp(i));
	
		g.vm_add_model("Sph", 1, 4);
		g.vm_add_model("Nug", 1, 0);
		g.vm_print();
	
		g.v_fit();
		printf("fitted: ");
		g.vm_print();
	}

	catch (GStatFatalError) {
		printf("caught a GstatFatalError!\n");
	}
	catch (GStatNormalError) {
		printf("caught a GstatNormalError!\n");
	}
	catch (...) {
		printf("caught an (...) error!\n");
	}

}
#endif
