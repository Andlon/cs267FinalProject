#include <qwt_plot.h>
#include "gui.h"

main(int argc, char *argv[]) {
	QwtPlot *myPlot;
	long curve1, curve2;                    // keys
	double x[100], y1[100], y2[100];        // x and y values
	
	myPlot = new QwtPlot("Two Graphs", parent, "variogram");
	
	// add curves
	curve1 = myPlot->insertCurve("Graph 1");
	curve2 = myPlot->insertCurve("Graph 2");
	
	getSomeValues(x, y1, y2);
	
	// copy the data into the curves
	myPlot->setCurveData(curve1, x, y1, 100);      
	myPlot->setCurveData(curve2, x, y2, 100);
	
	// finally, refresh the plot
	myPlot->replot();              
}
