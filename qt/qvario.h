#include <stdio.h>

#include <qapplication.h>
#include <qtooltip.h>
#include <qlineedit.h>
#include <qlayout.h>
#include <qpushbutton.h>
#include <qmessagebox.h>
#include <qslider.h>
#include <qlabel.h>

#include "data.h"
#include "vario.h"
#include "glvars.h"
#include "utils.h"

/*
 * class for modifying a a variogram model, VARIOGRAM (see vario.h)
 */
class QVario: public QWidget
{
	Q_OBJECT
private:
	QLineEdit *le;
	QSlider *nugget_sb, *sill_sb, *range_sb, *relnug_sb;
	VARIOGRAM *vptr, *v_tmp;
	int nugget_perc, sill_perc, range_perc, relative_nugget_perc,
		text_changed;
	double range_guess, sill_guess;

public: 
	QVario(VARIOGRAM *v);
	~QVario(void);	

private:
	void UpdateVariogram(void);
	void reset_perc(void);
	double guess_sill(void);
	double guess_range(void);
	void set_default_model(void);
	void get_percentages(int reset);
	
private slots:
	void OKClicked(void);
	void ApplyClicked(void);
	void CancelClicked(void);
	void NuggetSliderChanged(int);
	void SillSliderChanged(int);
	void RangeSliderChanged(int);
	void RelativeNuggetSliderChanged(int);
	void LineEditTextChanged(const char *str);

signals:
	void VariogramChanged(void);
};
