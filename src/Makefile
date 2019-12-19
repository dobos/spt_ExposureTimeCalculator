.c.o:
	$(CC) -c -o $*.o $(CFLAGS) $*.c

OPT := -O3
CFLAGS = $(OPT) -DHGCDTE_SUTR -DMOONLIGHT_
LDFLAGS = -lm
BIN = ../python/pfsspecsim/bin

gsetc.x : libgsetc.o gsetc.o
	$(CC) gsetc.o libgsetc.o $(LDFLAGS) -o $(BIN)/gsetc.x
gsetc.o :
	$(CC) gsetc.c $(LDFLAGS) -c -o gsetc.o
libgsetc.a :
	$(CC) libgsetc.c $(LDFLAGS) -c -o libgsetc.o

clean :
	$(RM) $(BIN)/gsetc.x *~
	$(RM) gsetc.o *~
	$(RM) libgsetc.o *~
