PART1 =	abstract introduction linear design what implementation \
		technical proj discussion credits references 

PART2 =	decis1.eps decis1.fig decis2.eps decis2.fig figs makefile \
		mydefs report simulation

ALL = $(PART1)

BIN = ../src/gstat -s

all:
	lout -a -o all.ps all.lt 

asc:
	lout -a -P -o all.txt all.lt 

spell:;	ispell -p ../dict $(PART1)

wc:;	wc $(PART1)

clean:; rm -f lout.li *\.l[di] *\.l[di]x *\.bak louti*

vi:;	vi $(ALL)

figs:;	lout -s -o figs.ps figs.lt

aanb:;	lout -a -s -o aanb.ps aanb.lt

aanb2:;	lout -a -s -o aanb2.ps aanb2.lt

check:;	$(BIN) frag1
		$(BIN) frag1a
		$(BIN) frag2
		$(BIN) -c frag3
		$(BIN) -c frag4
		$(BIN) -c frag5
		$(BIN) -c frag6
		grep -v variogram frag6 | $(BIN) -c -
		cat frag6 frag7 | $(BIN) -c -
		cat frag5 frag7 frag8 | $(BIN) -c -
		cat frag5 frag7 frag8 frag9 | $(BIN) -c -
