#include <qapp.h>
#include <qwidget.h>
#include <qprinter.h>

#include <qwt_plot.h>
#include <qwt_curve.h>

#include "qvario.h"
#include "gstatwrp.h"

class MainWin : public QWidget
{
	Q_OBJECT

private:
    QwtPlot *plot;
    Gstat *g;
    QVario *vario_modelw;

    long est, model, sel;
    double maxX, maxY, minY, selXmin, selXmax, selYmin, selYmax;
    int selecting;
    QPoint p1, p2;

public:
    MainWin(int argc, char *argv[]);
    ~MainWin();

protected:
    void resizeEvent(QResizeEvent *e);
	void keyPressEvent(QKeyEvent *k);
	void plot_model(long key);
	void plot_est(long key, long sel_key);

private slots:
	void plotMousePressed(const QMouseEvent &e);
	void plotMouseReleased(const QMouseEvent &e);
	void plotMouseMoved(const QMouseEvent &e);
	void legendClicked(long key);
	void plot_model(void);
};
