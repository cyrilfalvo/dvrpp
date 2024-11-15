#
#  Makefile 
#     make all   (construction de l'executable)
#     make clean (effacement des fichiers objets et de l'executable)
#

CC      = g++-14
CFLAGS   = -O2 -Wall 
CLIBS     = ../libdvr.a -llapack

FF      = gfortran-14
FFLAGS   = -cpp -O2 -Wall 
FLIBS     = ../libdvr.a -llapack -lstdc++

OBJ_DIR = obj
SRC_DIR = src

PROG_CPP     = ex_morse_cpp.x
PROG_FOR     = ex_morse_fortran.x

LISTEOBJ_CPP = \
	${OBJ_DIR}/ex_morse_cpp.o

LISTEOBJ_FOR = \
	${OBJ_DIR}/DVRpp_fortran_interface.o\
	${OBJ_DIR}/ex_morse_fortran.o


#
# compilation
#
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FF) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@
#
# dependances
#
all       : $(PROG_CPP) $(PROG_FOR)
#
# edition de lien
#
$(PROG_CPP) : $(LISTEOBJ_CPP)
	$(CC) -o $(PROG_CPP) $(CFLAGS) $(LISTEOBJ_CPP) $(CLIBS)

$(PROG_FOR) : $(LISTEOBJ_FOR)
	$(FF) -o $(PROG_FOR) $(FFLAGS) $(LISTEOBJ_FOR) $(FLIBS)
#
# effacement des fichiers intermediaires
#
clean :
	rm -f *.mod $(PROG_CPP) $(PROG_FOR) $(LISTEOBJ_CPP) $(LISTEOBJ_FOR)

