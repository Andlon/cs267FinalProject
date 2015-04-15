#ifndef GSTATWRP_H /* avoid multiple inclusion: */
# define GSTATWRP_H
/*! \file gstatwrp.h
    \brief gstat C++ wrapper functions
*/

extern "C" {
	extern void gstat_main(int *argc_p, char **argv_p[]);
}

#ifndef DATA_H
# define DATA void
#endif

#ifndef VARIO_H
# define VARIOGRAM void
#endif

#ifndef INCLUDED_APP_GSTAT
#include "gstatio.h"
#endif

/*! gstat class, initializes and activetes gstat. 

	Only one class gstat can be active at a time
*/
class Gstat {
	// private members:
	static int active; //!< prevent multiple simultaneous instances of Gstat
	DATA *dptr; //!< pointer to current data variable
	VARIOGRAM *vptr; //!< pointer to current variogram structure
	void add_data(int id); //!< read a parsed data structure 
	int add_edit_model(int id, const char *name, double sill, double range);
		//!< enter parameters for model name, add (id < 0) or edit (id >= 0)
	void process_gstat_error(int level);
		//!< do something with a gstat error condition

public:
	Gstat(app_GStat *pt); //!< constructor
	~Gstat(); //!< destructor
	const char *get_version() const; //!< returns version info: "2.0x (date)"
	void add_data(const char *name, const char *file_name,
			int x = 0, int y = 0, int v = 0); 
			//!< add a data variabble
	void remove_data(const char *name); //!< remove variable name 
	size_t get_nvars() const; //!< returns Gstat.nvars
	void report(FILE *f, const char *name) const; 
		//!< write report of var name to f
	void report(FILE *f, int id) const; 
		//!< write report of var id to f
	int read_commands(const char *cmd, const char *file);
		//!< parse gstat commands, read data variables (when present)
	void write_commands(const char *fname); 
		//!< write command file to fname
	const char *d_name(int id) const; 
		//!< get name of variable id
	size_t d_id(const char *name) const;
		//!< get id of variable name

	// d_* functions refer to data query functions
	void d_set(const char *name); 
		//!< set current data variable to name (only relevant when nvars > 1)
	void d_set(int id); 
		//!< set current data variable to id
	size_t d_get_n() const; 
		//!< number of data points in this variable
	void d_get_xyv(int i, double *x, double *y, double *v) const;
		//!< get x,y coordinates and value of data list number i

	// v_* functions refer to variogram query functions
	void v_set(const char *name); 
		//!< make variogram for variable name active
	void v_set(const char *a, const char *b); 
		//!< make cross variogram (a,b) active
	void v_calc(const char *name); 
		//!< calculate variogram for variable name
	void v_calc(const char *a, const char *b); 
		//!< calculate cross variogram for variables a and b 
	double v_get_model(double dx, double dy, double dz) const; 
		//!< get variogram model value for vector [dx,dy,dz]
	void v_fit(); 
		//!< fit current variogram model, if present
	size_t v_get_n_est() const;
		//!< get number of sample variogram estimates
	double v_get_est(size_t i) const;
		//!< get sample variogram value i
	double v_get_dist(size_t i) const;
		//!< get average distance for value i
	size_t v_get_npp(size_t i) const;
		//!< get number of pairs of points for value i
	void v_get_ppair(size_t i, size_t j,
			double *x1, double *y1, double *x2, double *y2) const; 
		//!< for j = 0 ... v_get_npp(i), return point pair (x1,y1),(x2,y2)
	void v_set_cutoff(double val);
		//!< set cutoff for semivariogram calculation
	double v_get_cutoff() const;
		//!< get cutoff for semivariogram calculation
	void v_set_width(double val); 
		//!< set interval width for semivariogram calculation
	double v_get_width() const; 
		//!< get interval width for semivariogram calculation

	// vm_* function to query/set variogram model stuff
	int    vm_has_model() const;
		//!< non-zero if a model is present
	int    vm_add_model(const char *name, double sill, double range);
		//!< add the named model, and return a unique id for this model
	void   vm_edit_model(int id, const char *name, double sill, double range);
		//!< edit (modify) the submodel with this id
	void   vm_delete_model(int id);
		//!< delete the submodel with this id
	const char **vm_get_possible_models() const; 
		//!< NULL-terminated char* array with all permitted models
	const char *vm_get_nugget_string() const;
		//!< return "Nug"
	const char *vm_get_model_description() const;
		//!< return a string with the model, e.g. `1.5 Nug() + 2.5 Exp(3.2)'
	double vm_get_default_range(const char *model) const;
		//!< return a ``suitable'' default range for this model
	double vm_get_default_sill() const;
		//!< return a ``suitable'' default sill 
	double vm_get_default_nugget() const; 
		//!< get a default sill value
	double vm_get_sill(int id) const;
		//!< get a sill value for model part id
	double vm_get_range(int id) const;
		//!< get range value of part id
	void   vm_print() const;
		//!< debug purposes: print model to stdout

	void set_warning_handler(void (*warning_fn)(const char *message));
		//!< set up the gstat warning handler (do not exit after a warning)
	void set_error_handler(void (*error_fn)(const char *message, int level));
		//!< set up the gstat error handler (should exit after an error)
};

#endif /* GSTATWRP_H */
