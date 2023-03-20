DSRC2=../src

F_COMP=gfortran
F_OPTS= -O0 -g -fcheck=bounds -fbacktrace -ffree-line-length-200 -Wall -fdefault-real-8

LOADER=$(F_COMP)
LOADER_OPTS=-v #-Bstatic -check_bounds -traceback #-Mlarge_arrays
ARCHIVE=ar rs

# Compiler commands.

F_COMMAND = $(F_COMP) -c $(F_OPTS) $(INCLUDES)

# Define archive and executable names.

BASE=gf
EXE=$(BASE).x
ARC=$(BASE).a

# Define source.
F_SRC=  $(DSRC2)/modConstants.F90 \
	$(DSRC2)/modHenrysLawCts.F90 \
	$(DSRC2)/modGate.F90 \
	$(DSRC2)/utils/modVector/modVector.F90 \
	$(DSRC2)/modConvParGF.F90 \
	$(DSRC2)/1d-gf-test.f90 	

# Define targets.

all: $(EXE)

$(EXE): $(ARC) 
	$(LOADER) -o $(EXE) $(LOADER_OPTS) $(ARC) $(LIB_IFS)

	rm -f r*.o

$(ARC): $(F_SRC)
	$(F_COMMAND) $(?)
	$(ARCHIVE) $(ARC) *.o

install:
	ln -fs $(pwd)/$(EXE) ../$(BASE)

clean:
	rm -f $(ARC) $(EXE) $(BASE) *.o  ___* core* *.mod *.a

modConstants.o : $(DSRC2)/modConstants.F90

modGate.o : $(DSRC2)/modGate.F90 modConstants.o

modHenrysLawCts.o : $(DSRC2)/modHenrysLawCts.F90 modConstants.o

modVector.o : $(DSRC2)/utils/modVector/modVector.F90 

modConvParGF.o : $(DSRC2)/modConvParGF.F90 modConstants.o modGate.o modHenrysLawCts.o modVector.o

1d-gf-test.o : $(DSRC2)/1d-gf-test.f90 modConstants.o modConvParGF.o


