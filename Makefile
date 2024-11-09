#
#  Makefile 
#     make all   (construction de l'executable)
#     make clean (effacement des fichiers objets et de l'executable)
#
 

CC       = icpc 
AR       = ar
CFLAGS   = -O3 -D_ARPACK_

OBJ_DIR = obj
SRC_DIR = src
#CC       = pgCC
#LD       = pgCC
#CFLAGS    = -tp p7-64 -fastsse -O4  -pgf90libs 
#LIBS      = -L${LIBDIR} -L/home/cyril_falvo/opt/ARPACK -larpack_LINUX -lacml -lfftw3
#

PROG     = libdvr.a
LISTEOBJ = \
        ${OBJ_DIR}/BasisSet.o \
        ${OBJ_DIR}/BasisHerm.o \
        ${OBJ_DIR}/BasisSine.o \
        ${OBJ_DIR}/BasisExp.o \
        ${OBJ_DIR}/DVR++.o \
#
# compilation
#
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@
#
# dependances
#
all       : $(PROG)
#
# edition de lien
#
$(PROG) : $(LISTEOBJ)
	$(AR) cr $(PROG) $(LISTEOBJ) 
#
# effacement des fichiers intermediaires
#
clean :
	rm -f $(PROG) $(LISTEOBJ)

