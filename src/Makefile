CC = gcc

CFLAGS = -O -g -Wall
LFLAGS = -lm

HEADER = constraints.h reader_energy_params.h util.h energy_functions.h sampling.h mfe_backtrack.h RNAmutants.h hybrid.h energy_params_tables.h
OBJ = constraints.o reader_energy_params.o util.o energy_functions.o sampling.o mfe_backtrack.o RNAmutants.o hybrid.o
EXE = RNAmutants

.SUFFIXES:
.SUFFIXES: .c .h .o .l $(EXE)

all :	RNAmutants.c RNAmutants

%.o :	%.c %.h
	$(CC) $(CFLAGS) -c $*.c

%$(EXE) : %.o
	$(CC) $(CFLAGS) $(LFLAGS) $? -o $@

RNAmutants : $(OBJ)
	$(CC) $(LFLAGS) $(CFLAGS) $(OBJ) -o $@

clean:	
	rm -f *~ $(OBJ)
