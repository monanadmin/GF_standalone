# Novo Makefile, preparado para compilar com gfortran e intel:
#
# Usage: make usage, ou apenas make
#
#			make gnu		:: para compilar com gfortran
#			make intel  :: para compilar com ifort
#
# Carlos Renato
# 2023/mar
#-------------------------------------------------------------------

BIN=../bin
SRC=../src
UTILS=$(SRC)/utils

BASE=gf
EXE=$(BASE).x
ARC=$(BASE).a

# Opcoes para o compilador gfortran:
GNU_OPTS=-O0 -g -fcheck=bounds -fbacktrace -ffree-line-length-200 -fdefault-real-8

# Opcoes para o compilador ifort:
INTEL_OPTS=-O0 -g -check bounds -traceback -check uninit -r8 -free

dummy:
	@$(MAKE) usage --no-print-directory

gnu:
	@( $(MAKE) all --no-print-directory \
	"COMP = gfortran" \
	"F_OPTS = $(GNU_OPTS)" )

intel:
	@( $(MAKE) all --no-print-directory \
	"COMP = ifort" \
	"F_OPTS = $(INTEL_OPTS)" )

all:	$(BIN)/$(EXE)

$(BIN)/$(EXE):	$(BIN)/$(ARC)
	$(COMP) $(F_OPTS) $< -o $@

$(BIN)/$(ARC):	$(BIN)/1d-gf-test.o $(BIN)/modConvParGF.o $(BIN)/modConstants.o $(BIN)/modGate.o $(BIN)/modHenrysLawCts.o $(BIN)/modVector.o 
	ar rs $@ $(BIN)/*.o

$(BIN)/1d-gf-test.o: $(SRC)/1d-gf-test.f90 $(BIN)/modGate.o $(BIN)/modConvParGF.o 
	$(COMP) $(F_OPTS) -c $< -o $@

$(BIN)/modConvParGF.o: $(SRC)/modConvParGF.F90 $(BIN)/modConstants.o $(BIN)/modGate.o $(BIN)/modHenrysLawCts.o $(BIN)/modVector.o
	$(COMP) $(F_OPTS) -c $< -o $@

$(BIN)/modConstants.o:	$(SRC)/modConstants.F90
	$(COMP) $(F_OPTS) -c $< -o $@

$(BIN)/modGate.o:	$(SRC)/modGate.F90
	$(COMP) $(F_OPTS) -c $< -o $@

$(BIN)/modHenrysLawCts.o:	$(SRC)/modHenrysLawCts.F90
	$(COMP) $(F_OPTS) -c $< -o $@

$(BIN)/modVector.o:	$(UTILS)/modVector/modVector.F90
	$(COMP) $(F_OPTS) -c $< -o $@

clean:
	rm -f $(BIN)/*.x $(BIN)/*.o $(BIN)/*.a $(BIN)/*.mod

usage:
	@echo ""
	@echo "Usage: $(MAKE) <target> "
	@echo ""
	@echo "Example targets:"
	@echo "                 gnu"
	@echo "                 intel"
	@echo ""
