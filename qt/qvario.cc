#include <assert.h>
#include <math.h>
#include "qvario.h"
#include "lex.h"

QVario::QVario(VARIOGRAM *v)
{
	assert(v != NULL);
	assert(v->id1 == v->id2);

	vptr = v; /* keep a local copy */
	if (vptr->n_models == 0)
		set_default_model();

	v_tmp = init_variogram(NULL); /* keep a scratch variogram */
	range_guess = 1000; // guess_range();
	sill_guess = guess_sill();
	relative_nugget_perc = (int) floor(100 * relative_nugget(vptr) + 0.5);
	get_percentages(0);

	// Create a layout to position the widgets
	QVBoxLayout *topLayout = new QVBoxLayout(this, 10);

	QGridLayout *grid = new QGridLayout(4, 3);
	// This layout will get all of the stretch
	topLayout->addLayout(grid); //, 10);

	// Create range slider
	QLabel *range_msgLabel = new QLabel(this, "msgLabel");
	range_msgLabel->setText("Range:");
	range_msgLabel->setAlignment(AlignRight|AlignVCenter);
	// range_msgLabel->setFixedSize(msgLabel->sizeHint());
	grid->addWidget(range_msgLabel, 0, 0, AlignRight|AlignVCenter);
	range_sb = new QSlider(1, 300, 1, range_perc, QSlider::Horizontal,
			this, "Slider");
	range_sb->setTickmarks(QSlider::Below);
	range_sb->setTickInterval(10);
	range_sb->setFocusPolicy(QWidget::TabFocus);
	range_sb->setFixedHeight(range_sb->sizeHint().height());
	// vbox->addWidget(range_sb);
	grid->addMultiCellWidget(range_sb, 0, 0, 1, 1);
	connect(range_sb, SIGNAL(valueChanged(int)), SLOT(RangeSliderChanged(int)));
	// QToolTip::add(range_sb, "adjust range parameter");

	// Create nugget slider
	QLabel *nugget_msgLabel = new QLabel(this, "msgLabel");
	nugget_msgLabel->setText("Nugget:");
	nugget_msgLabel->setAlignment(AlignRight|AlignVCenter);
	// nugget_msgLabel->setFixedSize(msgLabel->sizeHint());
	grid->addWidget(nugget_msgLabel, 1, 0, AlignVCenter);
	nugget_sb = new QSlider(1, 200, 1, nugget_perc,
			 QSlider::Horizontal, this, "Slider");
	nugget_sb->setTickmarks(QSlider::Below);
	nugget_sb->setTickInterval(10);
	nugget_sb->setFocusPolicy(QWidget::TabFocus);
	nugget_sb->setFixedHeight(nugget_sb->sizeHint().height());
	// vbox->addWidget(sb);
	grid->addMultiCellWidget(nugget_sb, 1, 1, 1, 1);
	connect(nugget_sb, SIGNAL(valueChanged(int)),
			SLOT(NuggetSliderChanged(int)));
    // QToolTip::add(nugget_sb, "adjust nugget");

	QLabel *relnug_msgLabel = new QLabel(this, "msgLabel");
	relnug_msgLabel->setText("Relative Nugget:");
	relnug_msgLabel->setAlignment(AlignRight|AlignVCenter);
	// relnug_msgLabel->setFixedSize(msgLabel->sizeHint());
	grid->addWidget(relnug_msgLabel, 2, 0, AlignVCenter);

	relnug_sb = new QSlider(0, 99, 1, relative_nugget_perc, QSlider::Horizontal,
			this, "Slider");
	relnug_sb->setTickmarks(QSlider::Below);
	relnug_sb->setTickInterval(3);
	relnug_sb->setFocusPolicy(QWidget::TabFocus);
	relnug_sb->setFixedHeight(relnug_sb->sizeHint().height());
	// vbox->addWidget(sb);
	grid->addMultiCellWidget(relnug_sb, 2, 2, 1, 1);
	connect(relnug_sb, SIGNAL(valueChanged(int)), 
		SLOT(RelativeNuggetSliderChanged(int)));
	// QToolTip::add(relnug_sb, "adjust relative nugget");

	// Create sill slider
	QLabel *sill_msgLabel = new QLabel(this, "msgLabel");
	sill_msgLabel->setText("Sill:");
	sill_msgLabel->setAlignment(AlignRight|AlignVCenter);
	// sill_msgLabel->setFixedSize(msgLabel->sizeHint());
	grid->addWidget(sill_msgLabel, 3, 0, AlignVCenter);
	sill_sb = new QSlider(1, 200, 1, sill_perc, QSlider::Horizontal,
			this, "Slider");
	sill_sb->setTickmarks(QSlider::Below);
	sill_sb->setTickInterval(10);
	sill_sb->setFocusPolicy(QWidget::TabFocus);
	sill_sb->setFixedHeight(sill_sb->sizeHint().height());
	// vbox->addWidget(sb);
	grid->addMultiCellWidget(sill_sb, 3, 3, 1, 1);
	connect(sill_sb, SIGNAL(valueChanged(int)), SLOT(SillSliderChanged(int)));
	// QToolTip::add(sill_sb, "adjust sill");

	// reset_perc();

  // Create a horizontal line (sort of QFrame) above the message line

	QGridLayout *grd2 = new QGridLayout(3, 3);
	topLayout->addLayout(grd2, 0);

    QFrame *separator = new QFrame(this, "separatorLine");
    separator->setFrameStyle(QFrame::HLine | QFrame::Sunken);
    separator->setFixedHeight(separator->sizeHint().height());
    grd2->addMultiCellWidget(separator, 0, 0, 0, 2);

	// Create a single line edit
	le = new QLineEdit(this, "lineEdit");
	grd2->addMultiCellWidget(le, 1, 1, 0, 2);
	// le->setFixedHeight(le->sizeHint().height());
	le->setFixedHeight(20);
	le->setText(v->descr);
	connect(le, SIGNAL(textChanged(const char *)), 
			SLOT(LineEditTextChanged(const char *)));
	connect(le, SIGNAL(returnPressed()), SLOT(ApplyClicked()));
	// QToolTip::add(le, "Edit text to adjust variogram model");
	text_changed = 0;

	// Create Apply Push button
	QPushButton *pb_apply;
	pb_apply = new QPushButton(this, "Apply button");    // create button 1
	pb_apply->setText("Apply");
	// pb_apply->setFixedHeight(pb_apply->sizeHint().height());
	pb_apply->setMinimumSize(pb_apply->sizeHint());
	grd2->addWidget(pb_apply, 2, 0, AlignVCenter);
	// QToolTip::add(pb_apply, "Apply changes to plot");
	connect(pb_apply, SIGNAL(clicked()), SLOT(ApplyClicked()));

	// Create OK Push button
	QPushButton *pb_ok;
	pb_ok = new QPushButton(this, "OK button");    // create button 1
	pb_ok->setText("OK");
	// pb_ok->setFixedHeight(pb_ok->sizeHint().height());
	pb_ok->setMinimumSize(pb_ok->sizeHint());
	grd2->addWidget(pb_ok, 2, 1, AlignVCenter);
	// QToolTip::add(pb_ok, "Apply changes to plot and close");
	connect(pb_ok, SIGNAL(clicked()), SLOT(OKClicked()));

	// Create Cancel Push button
	QPushButton *pb_cancel;
	pb_cancel = new QPushButton(this, "Cancel button");    // create button 1
	pb_cancel->setText("Cancel");
	// pb_cancel->setFixedHeight(pb_cancel->sizeHint().height());
	pb_cancel->setMinimumSize(pb_cancel->sizeHint());
	grd2->addWidget(pb_cancel, 2, 2, AlignVCenter);
	// QToolTip::add(pb_cancel, "Abandon changes and close");
	connect(pb_cancel, SIGNAL(clicked()), SLOT(CancelClicked()));

	grid->setRowStretch(0,1);
	grid->setRowStretch(1,1);
	grid->setRowStretch(2,1);
	grid->setRowStretch(3,1);

	grid->setColStretch(0,1);
	grid->setColStretch(1,4);
	grid->setColStretch(2,1);

	// grd2->setColStretch(0,0);
	// grd2->setColStretch(1,0);
	// grd2->setColStretch(2,0);

	// topLayout->activate();
	topLayout->addStrut(0);
}

QVario::~QVario(void)
{
	free_variogram(v_tmp);
	return;
}

void QVario::LineEditTextChanged(const char *str)
{
	if (str != NULL)
		text_changed = 1;
}

void QVario::OKClicked(void)
{
	UpdateVariogram();
	this.hide();
}

void QVario::CancelClicked(void)
{
	le->setText(vptr->descr);
	text_changed = 0;
#ifdef QVARIO_MAIN
	QApplication::exit(0);
#else
	this.hide();
#endif
}

void QVario::ApplyClicked(void)
{
	text_changed = 1;
	UpdateVariogram();
}

void QVario::RangeSliderChanged(int val)
{
	if (val != range_perc) 
	{
		double dval = 1.0 * val;
		for (int i = 0; i < vptr->n_models; i++)
			vptr->part[i].range *= dval/range_perc;
		range_perc = val;
		update_variogram(vptr);
		le->setText(vptr->descr);
		text_changed = 0;
		emit VariogramChanged();
	}
}

void QVario::SillSliderChanged(int val)
{
	int i;

	if (val != sill_perc) 
	{
		// val is the new percentage, relative to sill_guess
		// modify only non-nugget sills, proportionally, to get there.
		double sum_sill = 0.0, sum_nug = 0.0, fraction = 0.0;
		double dval = val / 100.0;
		for (i = 0; i < vptr->n_models; i++)
			if (vptr->part[i].model != NUGGET)
				sum_sill += vptr->part[i].sill;
			else
				sum_nug += vptr->part[i].sill;
		assert(sum_sill > 0.0);
		if (dval * sill_guess > 1.01 * sum_nug)
		{
			fraction = (dval * sill_guess - sum_nug)/(sum_sill);

			for (i = 0; i < vptr->n_models; i++)
				if (vptr->part[i].model != NUGGET)
					vptr->part[i].sill *= fraction;

			sill_perc = (int) floor(dval * 100 + 0.5);
		}

		update_variogram(vptr);
		le->setText(vptr->descr);
		text_changed = 0;
		emit VariogramChanged();
	}
}

void QVario::NuggetSliderChanged(int val)
{
	if (val != nugget_perc) 
	{
		double dval = 1.0 * val;
		for (int i = 0; i < vptr->n_models; i++)
			if (vptr->part[i].model == NUGGET)
				vptr->part[i].sill = dval * sill_guess/nugget_perc;
		nugget_perc = val;

		update_variogram(vptr);
		le->setText(vptr->descr);
		text_changed = 0;
		emit VariogramChanged();
	}
}

void QVario::RelativeNuggetSliderChanged(int val)
{
	if (vptr->n_models > 1)
	{
		double dval = (1.0 * val)/100.0;
		double nug = 0.0, sum = 0.0;
		for (int i = 0; i < vptr->n_models; i++)
		{
			sum += vptr->part[i].sill;
			if (vptr->part[i].model == NUGGET)
				nug += vptr->part[i].sill;
		}
	
		for (int i = 0; i < vptr->n_models; i++)
		{
			if (vptr->part[i].model == NUGGET)
				vptr->part[i].sill = dval * sum;
			else
				vptr->part[i].sill =
					(1.0 - dval) * sum * (vptr->part[i].sill/(sum-nug));
		}
		relative_nugget_perc = val;
		update_variogram(vptr);
		le->setText(vptr->descr);
		text_changed = 0;
		emit VariogramChanged();
	}
}

void QVario::reset_perc(void)
{
	nugget_sb->setValue(nugget_perc = 100);
	sill_sb->setValue(sill_perc = 100);
	range_sb->setValue(range_perc = 100);
	relative_nugget_perc = (int) floor(100 * relative_nugget(vptr) + 0.5);
	relnug_sb->setValue(relative_nugget_perc);
}

void QVario::get_percentages(int reset)
{
	int i;
	double sum_s = 0.0, sum_n = 0.0, max_range = 0.0;

	for (i = 0; i < vptr->n_models; i++)
	{
		if (vptr->part[i].model != NUGGET)
			sum_s += vptr->part[i].sill;
		else
			sum_n += vptr->part[i].sill;
		max_range = MAX(max_range, vptr->part[i].range);
	}
	sill_perc = (int) floor(100 * sum_s / sill_guess);
	nugget_perc = (int) floor(100 * sum_n / sill_guess);
	range_perc = (int) floor(100 * max_range / range_guess);

	if (reset) 
	{
		nugget_sb->setValue(nugget_perc);
		sill_sb->setValue(sill_perc);
		range_sb->setValue(range_perc);
	}

	printf("sill: %d, nugget: %d, range: %d\n", sill_perc, nugget_perc, range_perc);
}

void QVario::set_default_model(void)
{
	VGM_MODEL vm;
	vm.fit_sill = vm.fit_range = 0;

	// First structure:
	vm.sill = 0.3 * range_guess;
	vm.range = 0.0;
	vm.model = NUGGET;
	push_variogram_model(vptr, vm);

	// Second structure:
	vm.sill = 0.7 * range_guess;
	vm.range = range_guess;
	vm.model = SPHERICAL;
	push_variogram_model(vptr, vm);
}

double QVario::guess_range(void)
{
	if (vptr->id1 >= 0)
	{
		DATA **dpp = get_data();
		DATA *dp = dpp[vptr->id1];
		assert(dp);
		return 0.111 * data_block_diagonal(dp);
		// Why .111? --it's about one third of the default variogram distance
	}
	else
	{
		for (int i = 0; i < vptr->n_models; i++)
			if (vptr->part[i].range > 0.0)
				return vptr->part[i].range;
		return 1.0;
	}
}

double QVario::guess_sill(void)
{
	if (vptr->id1 >= 0)
	{
		assert(vptr->id1 == vptr->id2);
		DATA **d = get_data();
		int id = vptr->id1;
		assert(d[id]);
		return d[id]->std * d[id]->std;
	}
	else
	{
		if (vptr->n_models > 0)
			return vptr->sum_sills;
		else
			return 1.0;
	}
}

void QVario::UpdateVariogram(void)
{
	if (text_changed)
	{
		printf("reading `%s'\n", le->text());
		if (read_variogram(v_tmp, le->text())) {
			QMessageBox mb;
			mb.warning(this, "gstat",
				"Syntax error in variogram model,\nplease correct.");
		} else {
			read_variogram(vptr, le->text());
			update_variogram(vptr);
			vptr->ev->refit = 1;
			le->setText(vptr->descr);
			// reset_perc();
			get_percentages(1); // invokes a change due to roundoff --
			read_variogram(vptr, le->text());
			update_variogram(vptr);
			emit VariogramChanged();
		}
	}
	text_changed = 0;
}

#ifdef QVARIO_MAIN
int main(int argc, char **argv )
{
	VARIOGRAM *vp = init_variogram(NULL);
	if (argc > 1)
		read_variogram(vp, argv[1]);
	QApplication a(argc, argv);
	QVario w(vp);
	a.setMainWidget(&w);
	w.show();
	return a.exec();
}
#endif
