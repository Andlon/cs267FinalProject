C	= cgic.c cgictest.c ossfimcgi.c
H	= cgic.h
HTML = ossfim.html ossfim_form.html
DIV	= makefile
GIF = map.gif spacing.gif block.gif
WG = /home/edzer/c/whirlgif304/whirlgif

ossfimcgi:	cgic.o ossfimcgi.c
	gcc -c ossfimcgi.c
	gcc -o ossfimcgi ossfimcgi.o cgic.o

cgic.o:	cgic.c cgic.h
	gcc -c cgic.c

rcp:
	rcp -p $(C) $(H) $(HTML) $(DIV) gstat@gstat:c/cgi
	rcp -p $(HTML) $(GIF) gstat@gstat:data

install:	ossfimcgi
	cp ossfimcgi /www/gstat/cgi-bin

map.gif:	b4.gif
	cp b4.gif map.gif
	
b4.gif:	plt
	gnuplot plt

anim:
	$(WG) -time 75 -loop 0 s2.gif s3.gif s4.gif s5.gif s6.gif s7.gif \
		s8.gif s7.gif s6.gif s5.gif s4.gif s3.gif > spacing.gif
	$(WG) -time 75 -loop 0 b1.gif b2.gif b3.gif b4.gif b5.gif b6.gif \
		b5.gif b4.gif b3.gif b2.gif > block.gif
