#include <stdio.h>
#include <stdlib.h>

#define NLAGS 500

#include <ctype.h>
#include "utils.h" /* string_file() */
#include "data.h"
#include "vario.h"

#include "gstatwrp.h"
#include "gui_plot.h"

//-----------------------------------------------------------------
//		gui_plot.cc
//-----------------------------------------------------------------

MainWin::MainWin(int argc, char *argv[])
: QWidget(0,0)
{
    int i, n;

    if (argc <= 1) {
    	printf("usage: %s gstat_cmd_file\n", argv[0]);
    	exit(0);
    }

	g = new Gstat(argc, argv); // starts gstat
    g->read_commands(string_file(argv[1]), argv[1]);

    if (g->get_nvars() < 1) {
    	// do something useful, or:
    	printf("no variables read!\n");
    	exit(1);
    }

    // Allocate the QwtPlot widget 
    plot = new QwtPlot(g->d_name(0), this);

    if (! plot) {
    	printf("MainWin error!\n");
    	exit(-1);
    }

    //
    //  Calculate some values
    //
	g->v_set(g->d_name(0));
	g->v_calc(g->d_name(0));

	vario_modelw = new QVario(g->vptr);

    selXmin = selXmax = selYmin = selYmax = -9999.0;

    est = plot->insertCurve("sample variogram");
    QwtSymbol estSym;
	estSym.setStyle(QwtSymbol::Cross);
	estSym.setPen(black);
	estSym.setSize(5);
    plot->setCurveSymbol(est, estSym);
    plot->setCurveStyle(est, QwtCurve::NoCurve);

    sel = plot->insertCurve("selection");
	estSym.setPen(red);
   	plot->setCurveSymbol(sel, estSym);
    plot->setCurveStyle(sel, QwtCurve::NoCurve);

	plot_est(est, sel);

	model =	plot->insertCurve("variogram model");
	plot_model(model);

    // Insert markers
    
    //  ...a horizontal line at y = 0...
    // mY = plot->insertLineMarker("y = 0", QwtPlot::yLeft);     
    // plot->setMarkerYPos(mY, 0.0);

    //  ...a vertical line at x = 2 * pi
    // mX = plot->insertLineMarker("x = 2 pi", QwtPlot::xBottom);
    // plot->setMarkerXPos(mX, 6.284);

	//  plot->setMarkerYPos(plot->insertLineMarker("", QwtPlot::yLeft), 0.0);
    //	plot->setMarkerXPos(plot->insertLineMarker("", QwtPlot::xBottom), 0.0);
    
    // We want a legend
    // plot->setOutlineStyle(Qwt::Rect);
    // plot->enableOutline(TRUE);
    plot->enableLegend(TRUE);
    plot->enableGridX(FALSE);
    plot->enableGridY(FALSE);

    // Set axis titles
    plot->setAxisTitle(QwtPlot::xBottom, "distance");
    plot->setAxisTitle(QwtPlot::yLeft, "semivariance");
    plot->enableOutline(TRUE);

    // arrays can be deleted safely, since setCurveData 
    //  has copied the values.
    connect(plot, SIGNAL(plotMouseMoved(const QMouseEvent&)),
        SLOT(plotMouseMoved( const QMouseEvent&)));
    connect(plot, SIGNAL(plotMousePressed(const QMouseEvent &)),
        SLOT(plotMousePressed( const QMouseEvent&)));
    connect(plot, SIGNAL(plotMouseReleased(const QMouseEvent &)),
        SLOT(plotMouseReleased( const QMouseEvent&)));
	connect(plot, SIGNAL(legendClicked(long)), 
		SLOT(legendClicked(long)));
	connect(vario_modelw, SIGNAL(VariogramChanged(void)), 
		SLOT(plot_model(void)));
}

MainWin::~MainWin()
{
    delete plot;
    delete g;
}

void MainWin::plot_model(long key)
{
    int i;

    double *d = new double[NLAGS];
    double *z = new double[NLAGS];

    for (i = 0; i < NLAGS; i++) 
    {
    	d[i] = (i == 0 ? 1e-10 : i * maxX/NLAGS);
    	z[i] = g->v_get_model(d[i], 0.0, 0.0);
    }
    plot->setCurvePen(key, QPen(blue));
    plot->setCurveData(key, d, z, NLAGS);
    delete[] d;
    delete[] z;

    plot->replot();
}

void MainWin::plot_model(void)
{
	plot_model(model);
}

void MainWin::plot_est(long key, long sel_key)
{
    int i, n = g->v_get_n_est(), m;
    double x, y;
    double *x1 = new double[n];
    double *y1 = new double[n];
    double *x2 = new double[n];
    double *y2 = new double[n];

    maxX = maxY = minY = 0.0;
    for (i = m = n = 0; i < g->v_get_n_est(); i++) {
		if (g->v_get_npp(i) > 0) {
			x = g->v_get_dist(i);
			maxX = MAX(x, maxX);
			y = g->v_get_est(i);
			maxY = MAX(y, maxY);
			minY = MIN(y, minY);
			if (selXmax < 0.0 || (
					x < selXmin || x > selXmax ||
					y < selYmin || y > selYmax)) {
				x1[m] = x;
				y1[m] = y;
				m++;
			} else {
				x2[n] = x;
				y2[n] = y;
				n++;
			}
		}
    }

	plot->setAxisScale(QwtPlot::xBottom, 0.0, 1.03 * maxX, 0.0);
	plot->setAxisScale(QwtPlot::yLeft, 1.03 * minY, 1.03 * maxY, 0.0);

    plot->setCurveData(key, x1, y1, m);
	plot->setCurveData(sel_key, x2, y2, n);

    delete[] x1;
    delete[] y1;
    delete[] x2;
    delete[] y2;

    plot->replot();
}

void MainWin::resizeEvent(QResizeEvent *e)
{
    plot->resize(e->size());
    plot->move(0,0);
}

void MainWin::keyPressEvent(QKeyEvent *k)
{
    switch ( tolower(k->ascii()) ) {
        case 'p': {                              // reload
    			QBrush br(red);
    			QPen pn(yellow);
    	
    			QwtSymbol sym1;
    			sym1.setBrush(br);
    			sym1.setPen(pn);
    			sym1.setSize(11);
	
    			QPrinter p;
	
    			if (p.setup(0)) {
	    			this->plot->print(p, QwtFltrDim(200));
	    		}
	    	}
            break;
        case 'q':                               // quit
            QApplication::exit();
            break;
    }
}

void MainWin::plotMouseMoved(const QMouseEvent &e)
{

	// find neares parameter:

	//g->set_range(x);	
	//g->set_sill(y);	

	if (! selecting)
		plot_model(model);
}

void MainWin::plotMousePressed(const QMouseEvent &e)
{
    // store position
    p1 = e.pos();

    double x = plot->invTransform(QwtPlot::xBottom, e.pos().x());
    double y = plot->invTransform(QwtPlot::yLeft, e.pos().y());
    
    if (e.button() == RightButton)
    {
    	// printf("r button: x = %g y = %g\n", x, y);
    	plot->setOutlineStyle(Qwt::Rect); 
		plot->enableOutline(TRUE);
    	selecting = 1;
    }
    if (e.button() == LeftButton)
    {
    	// printf("l button: x = %g y = %g\n", x, y);
		plot->enableOutline(FALSE);
    	selecting = 0;
    }
    // update cursor pos display
    plotMouseMoved(e);
}

void MainWin::plotMouseReleased(const QMouseEvent &e)
{
    // int x1, x2, y1, y2;
    int lim;
    
    if (selecting) {
    	p2 = e.pos();
    	double x1 = plot->invTransform(QwtPlot::xBottom, p1.x());
    	double y1 = plot->invTransform(QwtPlot::yLeft, p1.y());
    	double x2 = plot->invTransform(QwtPlot::xBottom, p2.x());
    	double y2 = plot->invTransform(QwtPlot::yLeft, p2.y());
    	selXmin = MIN(x1,x2);
    	selYmin = MIN(y1,y2);
    	selXmax = MAX(x1,x2);
    	selYmax = MAX(y1,y2);
    	// printf("[%g,%g,%g,%g]\n",selXmin, selXmax, selYmin, selYmax);
		plot_est(est, sel);
    	selecting = 0;
    } else {
    	selXmin = selXmax = selYmin = selYmax = -9999.0;
		plot_est(est, sel);
    }
    // some shortcuts
    //int axl= QwtPlot::yLeft, axr = QwtPlot::yRight, axb= QwtPlot::xBottom;
    
    //if (d_zoom && (!d_zoomActive))
    //{
	//d_zoomActive = 1;
	
	// Don't invert any scales which aren't inverted
	//x1 = qwtMin(p1.x(), e.pos().x());
	//x2 = qwtMax(p1.x(), e.pos().x());
	//y1 = qwtMin(p1.y(), e.pos().y());
	//y2 = qwtMax(p1.y(), e.pos().y());
	
	// limit selected area to a minimum of 11x11 points
	//lim = 5 - (y2 - y1) / 2;
	//if (lim > 0)
	//{
	 //   y1 -= lim;
	 //   y2 += lim;
	//}
	//lim = 5 - (x2 - x1 + 1) / 2;
	//if (lim > 0)
	//{
	 //   x1 -= lim;
	 //  x2 += lim;
	//}
	
	// Set fixed scales
	// plt->setAxisScale(axl, plt->invTransform(axl,y1), plt->invTransform(axl,y2));
	// plt->setAxisScale(axr, plt->invTransform(axr,y1), plt->invTransform(axr,y2));
	// plt->setAxisScale(axb, plt->invTransform(axb,x1), plt->invTransform(axb,x2));
	// plt->replot();
	

    // }
    // frmInp->lblInfo->setText(cursorInfo);
    // plt->setOutlineStyle(Qwt::Triangle);
    
}

void MainWin::legendClicked(long key) 
{
	printf("key %d clicked\n", key);
	if (key == 3) {
		printf("here we go!\n");
		vario_modelw->show();
	}
}

int main(int argc, char **argv)
{

    QApplication a(argc, argv);
    MainWin w(argc, argv);
    a.setMainWidget(&w);
    w.show();
    return a.exec(); 
}
