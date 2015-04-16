CFLAGS = -Wall -g -I../gstatinc -I../meschach -Wall
LIBS = -L../lib -lgstat -lmeschach -lm -lcsf
# LIBS = -L /usr/local/qt/lib -L~/qt/qwt/lib -L/usr/X11R6/lib -L../lib \
#	 -lX11 -lgstat -lcsf -lmeschach -lm
CC = c++

gstatwrp:	../lib/libgstat.a gstatwrp.o gstatexcept.o
	$(CC) -o $@ $(CFLAGS) $(DEFS) gstatwrp.o $(LIBS)

gstatwrp.o:	gstatwrp.h gstatwrp.cc ../src/data.h ../src/glvars.h \
		../src/vario.h ../src/utils.h ../src/defs.h
	$(CC) -c $(CFLAGS) $(DEFS) gstatwrp.cc

gstatexcept.o:	gstatexcept.h gstatexcept.cc
	$(CC) -c $(CFLAGS) $(DEFS) gstatexcept.cc

gstatio.o:	gstatio.c gstatio.h
	gcc -c $(CFLAGS) $(DEFS) gstatio.c

errorwrp.o:	errorwrp.c errorwrp.h
	gcc -c $(CFLAGS) errorwrp.c

gui_plot:	$(OBJS) gstatwrp.o gui_plot.o gui_plot_moc.o qvario_moc.o qvario.o
	$(CC) -o $@ gui_plot.o gui_plot_moc.o gstatwrp.o qvario.o qvario_moc.o \
	$(LDFLAGS) $(LIBS) -lqt -lqwt

gui_plot.o:	gui_plot.cc gui_plot.h 
	$(CC) -c $(CFLAGS) -I/usr/local/qt/include \
		-I/home/edzer/qt/qwt/include gui_plot.cc

gui_plot_moc.o:	gui_plot.h
	moc -o gui_plot_moc.cc gui_plot.h
	$(CC) -c $(CFLAGS) -I/usr/local/qt/include -I/home/edzer/qt/qwt/include gui_plot_moc.cc

gstat-lib:	gstatwrp.o errorwrp.o gstatexcept.o gstatio.o
	ar ruv ../lib/libgstat.a gstatwrp.o errorwrp.o gstatexcept.o gstatio.o
	cp gstatwrp.h gstatexcept.h gstatio.h ../gstatinc

doc:
	doxygen doxygen.cfg

clean:
	rm -f gstatwrp gui_plot

qvario_moc.o:	qvario.h
	moc -o qvario_moc.cc qvario.h
	$(CC) -c $(CFLAGS) -I/usr/local/qt/include qvario_moc.cc

qvario.o:	qvario.h qvario.cc
	$(CC) -c $(CFLAGS) -I/usr/local/qt/include qvario.cc

qvario:	qvario_moc.o qvario.h qvario.cc
	$(CC) -c $(CFLAGS) -DQVARIO_MAIN -I/usr/local/qt/include qvario.cc
	$(CC) -o $@ qvario.o qvario_moc.o $(OBJS) $(LDFLAGS) $(LIBS) -lqt -lqwt
