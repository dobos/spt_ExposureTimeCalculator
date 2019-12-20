ifeq ($(CFG),debug)
	OPT = -g -Wall
else
	OPT = -O3
endif

CFLAGS = -DHGCDTE_SUTR -DMOONLIGHT_
LDFLAGS = -lm
BIN = ../python/pfsspecsim/bin
TARGETS = gssky gsetc.x

all : $(TARGETS)

gssky : libgsetc.o gssky.o
	$(CC) $(BIN)/gssky.o $(BIN)/libgsetc.o $(LDFLAGS) $(CFLAGS) $(OPT) -o $(BIN)/gssky
gssky.o :
	$(CC) gssky.c $(LDFLAGS) $(CFLAGS) $(OPT) -c -o $(BIN)/gssky.o	
gsetc.x : libgsetc.o gsetc.o
	$(CC) $(BIN)/gsetc.o $(BIN)/libgsetc.o $(LDFLAGS) $(CFLAGS) $(OPT) -o $(BIN)/gsetc.x
gsetc.o :
	$(CC) gsetc.c $(LDFLAGS) $(CFLAGS) $(OPT) -c -o $(BIN)/gsetc.o
libgsetc.o :
	$(CC) libgsetc.c $(LDFLAGS) $(CFLAGS) $(OPT) -c -o $(BIN)/libgsetc.o

clean :
	$(RM) $(BIN)/gsetc.x *~
	$(RM) $(BIN)/gsetc.o *~
	$(RM) $(BIN)/gssky *~
	$(RM) $(BIN)/gssky.o *~
	$(RM) $(BIN)/libgsetc.o *~
