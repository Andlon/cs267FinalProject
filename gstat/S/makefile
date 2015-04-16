VERSION	= 0.9-4

RSOURCE	= R/gstat_$(VERSION).tar.gz

SSOURCE	= gstat_Splus6_$(VERSION).tar.gz

RDEST	= R/gstat

SDEST	= Splus/gstat

SPLUS	= Splus61

CSRC	= 	../src/block.c \
	../src/block.h \
	../src/config.hin \
	../src/data.c \
	../src/data.h \
	../src/debug.h \
	../src/defaults.h \
	../src/defs.h \
	../src/direct.c \
	../src/direct.h \
	../src/fit.c \
	../src/fit.h \
	../src/getest.c \
	../src/getest.h \
	../src/gls.c \
	../src/gls.h \
	../src/glvars.c \
	../src/glvars.h \
	../src/lex.c \
	../src/lex.h \
	../src/lex.l \
	../src/lm.c \
	../src/lm.h \
	../src/lm_type.h \
	../src/makefile.in \
	../src/mapio.c \
	../src/mapio.h \
	../src/msim.c \
	../src/msim.h \
	../src/nsearch.c \
	../src/nsearch.h \
	../src/parse.h \
	../src/parse.y \
	../src/parse.c \
	../src/plot.c \
	../src/plot.h \
	../src/polygon.c \
	../src/polygon.h \
	../src/pqueue.c \
	../src/pqueue.h \
	../src/predict.c \
	../src/predict.h \
	../src/random.c \
	../src/random.h \
	../src/read.c \
	../src/read.h \
	../src/reml.c \
	../src/reml.h \
	../src/report.c \
	../src/report.h \
	../src/s.c \
	../src/select.c \
	../src/select.h \
	../src/sem.c \
	../src/sem.h \
	../src/sim.c \
	../src/sim.h \
	../src/stat.c \
	../src/stat.h \
	../src/userio.c \
	../src/userio.h \
	../src/utils.c \
	../src/utils.h \
	../src/vario.c \
	../src/vario.h \
	../src/vario_fn.c \
	../src/vario_fn.h \
	../src/vario_io.c \
	../src/vario_io.h \
	../src/version.h \
	../src/writecmd.c \
	../src/writecmd.h \
	../src/xvalid.c \
	../src/xvalid.h

MSRC	= ../meschach/chfactor.c \
	../meschach/copy.c \
	../meschach/err.c \
	../meschach/err.h \
	../meschach/hsehldr.c \
	../meschach/init.c \
	../meschach/ivecop.c \
	../meschach/lufactor.c \
	../meschach/machine.c \
	../meschach/machine.hin \
	../meschach/matop.c \
	../meschach/matrix.h \
	../meschach/matrix2.h \
	../meschach/matrixio.c \
	../meschach/meminfo.c \
	../meschach/meminfo.h \
	../meschach/memory.c \
	../meschach/memstat.c \
	../meschach/norm.c \
	../meschach/otherio.c \
	../meschach/pxop.c \
	../meschach/qrfactor.c \
	../meschach/solve.c \
	../meschach/sparse.c \
	../meschach/sparse.h \
	../meschach/sparse2.h \
	../meschach/sparseio.c \
	../meschach/spbkp.c \
	../meschach/spchfctr.c \
	../meschach/splufctr.c \
	../meschach/sprow.c \
	../meschach/spswap.c \
	../meschach/submat.c \
	../meschach/vecop.c \
	../meschach/version.c \
	../meschach/zmatrix.h \
	../meschach/zmatrix2.h

MMSRC	= ../meschach/dmacheps.c \
	../meschach/fmacheps.c \
	../meschach/maxint.c

default:
	make r-copy r-build
	R CMD INSTALL -l /home/edzer/S/library $(RSOURCE)

r-copy:
	cp -ra data demo man tests $(RDEST)
	cp -a $(CSRC) $(MSRC) $(RDEST)/src
	cp -a $(MMSRC) $(RDEST)/src/meschach
	./scripts/qtoR src $(RDEST)/R src/*q

r-build:
	(cd R; R CMD build gstat --force)

r-clean:
	rm -f $(RDEST)/src/*\.[chyl] 
	rm -f $(RDEST)/src/*in 
	rm -f $(RDEST)/src/meschach/* $(RDEST)/R/*
	rm -fr $(RDEST)/data $(RDEST)/demo $(RDEST)/man $(RDEST)/tests 
	rm -f $(RDEST)/*.tar.gz
	rm -fr R/gstat.Rcheck*
	(cd $(RDEST); ./cleanup)

r-check:
	(cd R; R CMD check gstat)

r-tgz:
	make r-clean r-copy r-build

r-rcp:	R/gstat.pdf R/gstat.dvi
	rcp $(RSOURCE) R/gstat.pdf R/gstat.dvi gstat@gstat.geog.uu.nl:data/bin

s-copy:
	cp -a $(CSRC) $(MSRC) $(SDEST)
	./scripts/preproc.pl -v SP5 ../src/s.c > $(SDEST)/s.c
	cp -a $(MMSRC) Splus/meschach
	./scripts/toSplus $(SPLUS)

s-build:
	(cd Splus; ./configure)
	(cd Splus; SPLUS=$(SPLUS) make library)

s-tgz:
	make s-clean s-copy
	(cd Splus; tar zcvf ../$(SSOURCE) *)

s-clean:
	(cd $(SDEST); rm -f *.[chyl] *.hin *.sgml *.q)
	rm -f Splus/config.status Splus/config.cache Splus/config.log 
	rm -f Splus/maxint Splus/macheps
	rm -fr $(SDEST)/.Data
	rm -fr $(SDEST)/*

s-rcp:
	make s-tgz
	rcp $(SSOURCE) gstat@gstat:data/bin

tgz:
	make r-tgz
	make s-tgz

rcp:	
	make r-rcp s-rcp

R/gstat.pdf:	./man
	rm -f R/gstat.pdf
	(cd R; R CMD Rd2dvi --pdf gstat)

R/gstat.dvi:	./man
	rm -f R/gstat.dvi
	(cd R; R CMD Rd2dvi gstat)

test.tgz:	test
	tar zcvf test.tgz test

flop:	test.tgz
	mcopy $(SOURCE) test.tgz makefile a:

clean:
	make r-clean s-clean

veryclean:
	make clean
	rm -f $(RSOURCE) $(SSOURCE) R/gstat.dvi R/gstat.pdf
