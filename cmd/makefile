#
# default target: process all non-interactive command files:
#
MAPS = sqrtdist part_a part_b mask_map

RESULT_MAPS = id_pr lzn_ckpr lzn_ckvr lzn_cspr lzn_okbp lzn_okbv lzn_okpr \
	lzn_okvr lzn_skpr lzn_skvr lzn_stpr lzn_stvr lzn_tr0 lzn_tr1 lzn_tr2 \
	lzn_tr3 lzn_trpr lzn_trvr lzn_ukpr lzn_ukvr lzn_uspr lzn_vkpr lzn_vkvr \
	lzn_vr0 lzn_vr1 lzn_vr2 lzn_vr3 lzn_stup lzn_stuv \
	sqd_ckpr sqd_ckvr sqd_okpr sqd_okvr znsqdcov

RESULT_ASCII = zinc_ok.out

BINDIR = ../src/

REFDIR = ../dos

BIN = gstat -s

GSTAT = $(BINDIR)$(BIN)

CMD = ex01.cmd ex02.cmd ex03.cmd ex04.cmd ex05.cmd ex06.cmd ex07.cmd \
	ex08.cmd ex09.cmd ex10.cmd ex11.cmd ex12.cmd ex13.cmd ex14.cmd ex14a.cmd \
	ex15.cmd

all:	$(MAPS)
		$(GSTAT) ex03.cmd
		$(GSTAT) ex04.cmd
		$(GSTAT) ex05.cmd
		$(GSTAT) ex06.cmd
		$(GSTAT) ex07.cmd
		$(GSTAT) ex08.cmd
		$(GSTAT) ex09.cmd
		$(GSTAT) ex10.cmd
		$(GSTAT) ex11.cmd
		$(GSTAT) ex12.cmd
		$(GSTAT) ex13.cmd
		$(GSTAT) ex14.cmd
		$(GSTAT) ex14a.cmd
		$(GSTAT) ex15.cmd
		$(GSTAT) ex16.cmd
		$(GSTAT) ex17.cmd

ex01:;	$(BINDIR)gstat ex01.cmd

ex02:;	$(BINDIR)gstat ex02.cmd

ex03:;	$(BINDIR)gstat ex03.cmd

ex04:;	$(BINDIR)gstat ex04.cmd

ex05:;	$(BINDIR)gstat ex05.cmd

ex06:;	$(BINDIR)gstat ex06.cmd

ex07:;	$(BINDIR)gstat ex07.cmd

ex08:;	$(BINDIR)gstat ex08.cmd

ex09:;	$(BINDIR)gstat ex09.cmd

ex10:;	$(BINDIR)gstat ex10.cmd

ex11:;	$(BINDIR)gstat ex11.cmd

ex12:;	$(BINDIR)gstat ex12.cmd

ex13:;	$(BINDIR)gstat ex13.cmd

ex14:;	$(BINDIR)gstat ex14.cmd

ex14a:;	$(BINDIR)gstat ex14a.cmd

ex15:;	$(BINDIR)gstat ex15.cmd

result_maps:
	grep predictions ex*cmd | awk -F "'" '{print $$2}'  > results
	grep variances ex*cmd | awk -F "'" '{print $$2}'  >> results
	grep covariances ex*cmd | awk -F "'" '{print $$2}'  >> results
	grep output ex*cmd | awk -F "'" '{print $$2}'  >> results

# remove files produced by test runs
clean:
	make result_maps
	rm -f `cat results` results $(RESULT_ASCII) fit.log

test:
	@echo "in case of equal maps, make test should result in Error 1"
	@sh map_diff $(REFDIR) $(RESULT_MAPS) | grep -v "are equal"

idrisi-ascii:
	rm -f $(MAPS)
	$(GSTAT) -e convert i part_a.csf part_a
	$(GSTAT) -e convert i part_b.csf part_b
	$(GSTAT) -e convert i mask_map.csf mask_map
	$(GSTAT) -e convert i sqrtdist.csf sqrtdist

idrisi-binary:
	rm -f $(MAPS)
	$(GSTAT) -e convert b part_a.csf part_a
	$(GSTAT) -e convert b part_b.csf part_b
	$(GSTAT) -e convert b mask_map.csf mask_map
	$(GSTAT) -e convert b sqrtdist.csf sqrtdist

arcfloat:
	rm -f $(MAPS)
	$(GSTAT) -e convert f part_a.csf part_a
	$(GSTAT) -e convert f part_b.csf part_b
	$(GSTAT) -e convert f mask_map.csf mask_map
	$(GSTAT) -e convert f sqrtdist.csf sqrtdist
	
arcgrid:
	rm -f $(MAPS)
	$(GSTAT) -e convert a part_a.csf part_a
	$(GSTAT) -e convert a part_b.csf part_b
	$(GSTAT) -e convert a mask_map.csf mask_map
	$(GSTAT) -e convert a sqrtdist.csf sqrtdist

maps.tgz:
	tar cf maps.tar $(RESULT_MAPS)
	gzip maps.tar
	mv maps.tar.gz maps.tgz
