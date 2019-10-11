FC = gfortran
FFLAGS = -ffast-math -Wall -march=native -O3 -mno-avx #compiler flags
LDFLAGS = -fopenmp #link flags

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

PROG = main #program name

#required objects: 
OBJS =
OBJS += constants.o
OBJS += plotroutines.o
OBJS += proc_output.o
OBJS += io.o
OBJS += initialize.o
OBJS += markov.o
OBJS += main.o

all: $(PROG)

main: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f95
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) *.mod
	$(RM) *.png *.txt *.plt *.dat
