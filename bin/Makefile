F90=gfortran

BIN=../bin
SRC=../src
UTILS=$(SRC)/utils


F_OPTS=-O0 -g -fcheck=bounds -fbacktrace -ffree-line-length-200 -fdefault-real-8 # -Wall -w  -std=f95 -> erro a partir do padrao 95

BASE=gf
EXE=$(BASE).x
ARC=$(BASE).a

all: $(BIN)/$(EXE)

$(BIN)/$(EXE):	$(BIN)/$(ARC)
	$(F90) $(F_OPTS) $< -o $@

$(BIN)/$(ARC):	$(BIN)/1d-gf-test.o $(BIN)/modConvParGF.o $(BIN)/modConstants.o $(BIN)/modGate.o $(BIN)/modHenrysLawCts.o $(BIN)/modVector.o 
	ar rs $@ $(BIN)/*.o

$(BIN)/1d-gf-test.o: $(SRC)/1d-gf-test.f90 $(BIN)/modGate.o $(BIN)/modConvParGF.o 
	$(F90) $(F_OPTS) -c $< -o $@

$(BIN)/modConvParGF.o: $(SRC)/modConvParGF.F90 $(BIN)/modConstants.o $(BIN)/modGate.o $(BIN)/modHenrysLawCts.o $(BIN)/modVector.o
	$(F90) $(F_OPTS) -c $< -o $@

$(BIN)/modConstants.o:	$(SRC)/modConstants.F90
	$(F90) $(F_OPTS) -c $< -o $@

$(BIN)/modGate.o:	$(SRC)/modGate.F90
	$(F90) $(F_OPTS) -c $< -o $@

$(BIN)/modHenrysLawCts.o:	$(SRC)/modHenrysLawCts.F90
	$(F90) $(F_OPTS) -c $< -o $@

$(BIN)/modVector.o:	$(UTILS)/modVector/modVector.F90
	$(F90) $(F_OPTS) -c $< -o $@


clean:
	rm -f $(BIN)/*.x $(BIN)/*.o $(BIN)/*.a $(BIN)/*.mod

	
