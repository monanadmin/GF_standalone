########################################################################
####################### Makefile Template ##############################
########################################################################

# Compiler settings - Can be customized.
CC = gfortran
DEBUG_FLAGS= -g -O0 -J obj -Wall -fcheck=all
FFLAGS?=${DEBUG_FLAGS}
LDFLAGS = 

# Makefile settings - Can be customized.
UNITY_TESTS?= Unity_tests
MODULAR_TESTS?= Modular_tests
UNITY_TESTS_APP_NAME?= Unity_tests_debug
MODULAR_TESTS_APP_NAME?= Modular_tests_debug
EXT = .F90
SRCDIR = src
TESTDIR = test
OBJDIR = obj

############## Do not change anything from here downwards! #############
SRC = $(wildcard $(SRCDIR)/*/*$(EXT))

OBJ = $(SRC:$(SRCDIR)/%$(EXT)=$(OBJDIR)/%.o)

MOD = *.mod
OBJTEST = $(TEST:$(TESTDIR)/%$(EXT)=$(TESTDIR)/%.o)

UNITY_TESTS_APP = ${OBJDIR}/test/${UNITY_TESTS_APP_NAME}
MODULAR_TESTS_APP = ${OBJDIR}/test/${MODULAR_TESTS_APP_NAME}

RM = rm -rf
########################################################################
####################### Targets beginning here #########################
########################################################################

all: $(UNITY_TESTS_APP) $(MODULAR_TESTS_APP)

# Builds the app Unity Tests
$(UNITY_TESTS_APP): ${TESTDIR}/${UNITY_TESTS}$(EXT) $(OBJ_TEST) $(OBJDIR)/model/moddata.o $(OBJDIR)/controller/modlist.o $(OBJDIR)/controller/modvector.o
	@mkdir -p $(OBJDIR)/test
	$(CC) $(FFLAGS) -o ${UNITY_TESTS_APP} $^ $(LDFLAGS)

# Builds the app Modular Tests
$(MODULAR_TESTS_APP): ${TESTDIR}/${MODULAR_TESTS}$(EXT)  $(OBJ_TEST) $(OBJDIR)/model/moddata.o $(OBJDIR)/controller/modlist.o $(OBJDIR)/controller/modvector.o
	$(CC) $(FFLAGS) -o ${MODULAR_TESTS_APP} $^ $(LDFLAGS)

$(OBJDIR)/model/moddata.o: $(SRCDIR)/model/moddata.F90
	@mkdir -p $(OBJDIR)/model
	$(CC) $(FFLAGS) -o $@ -c $<


$(OBJDIR)/controller/modlist.o: $(SRCDIR)/controller/modlist.F90 $(OBJDIR)/model/moddata.o
	@mkdir -p $(OBJDIR)/controller
	$(CC) $(FFLAGS) -o $@ -c $<


$(OBJDIR)/controller/modvector.o: $(SRCDIR)/controller/modvector.F90 $(OBJDIR)/model/moddata.o
	@mkdir -p $(OBJDIR)/controller
	$(CC) $(FFLAGS) -o $@ -c $<

	
################### Cleaning rules for Unix-based OS ###################
# Cleans complete project
.PHONY: clean

clean:
	$(RM) $(OBJ) $(OBJTEST) $(MOD) ${MODULAR_TESTS_APP} ${UNITY_TESTS_APP}

run: runUnityTest runModularTest

runUnityTest:
	(cd test && ./run_tests_with_type_of_param.sh ${UNITY_TESTS_APP_NAME})

runModularTest:
	(cd test && ./run_tests_with_type_of_param.sh ${MODULAR_TESTS_APP_NAME})

