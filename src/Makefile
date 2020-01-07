ifeq ($(CFG),debug)
	OPT = -g -Wall
else
	OPT = -O3
endif

CFLAGS = -DHGCDTE_SUTR -DMOONLIGHT_ -DNOPROGRESS
LDFLAGS = -lm
BIN = ../python/pfsspecsim/bin
TARGETS = gsnoise.x gspsf.x gsetc.x

all : $(TARGETS)

gsnoise.x : libgsetc.o gsnoise.o
	$(CC) $(BIN)/gsnoise.o $(BIN)/libgsetc.o $(LDFLAGS) $(CFLAGS) $(OPT) -o $(BIN)/gsnoise.x
gsnoise.o :
	$(CC) gsnoise.c $(LDFLAGS) $(CFLAGS) $(OPT) -c -o $(BIN)/gsnoise.o	
gspsf.x : libgsetc.o gspsf.o
	$(CC) $(BIN)/gspsf.o $(BIN)/libgsetc.o $(LDFLAGS) $(CFLAGS) $(OPT) -o $(BIN)/gspsf.x
gspsf.o :
	$(CC) gspsf.c $(LDFLAGS) $(CFLAGS) $(OPT) -c -o $(BIN)/gspsf.o	
gsetc.x : libgsetc.o gsetc.o
	$(CC) $(BIN)/gsetc.o $(BIN)/libgsetc.o $(LDFLAGS) $(CFLAGS) $(OPT) -o $(BIN)/gsetc.x
gsetc.o :
	$(CC) gsetc.c $(LDFLAGS) $(CFLAGS) $(OPT) -c -o $(BIN)/gsetc.o
libgsetc.o :
	$(CC) libgsetc.c $(LDFLAGS) $(CFLAGS) $(OPT) -c -o $(BIN)/libgsetc.o

clean :
	$(RM) $(BIN)/gsetc.x *~
	$(RM) $(BIN)/gsetc.o *~
	$(RM) $(BIN)/gspsf.x *~
	$(RM) $(BIN)/gspsf.o *~
	$(RM) $(BIN)/gsnoise.x *~
	$(RM) $(BIN)/gsnoise.o *~
	$(RM) $(BIN)/libgsetc.o *~
