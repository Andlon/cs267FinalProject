CC=g++
CFLAGS= -g -O3
DEPS = main.cpp
OBJ = main.o
 
default: serial

serial: $(OBJ)
	$(CC) $(CFLAGS) -o $@ $^

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp

clean:
	rm *.o serial 
