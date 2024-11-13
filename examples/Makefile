#
#  Makefile 
#     make all   (construction de l'executable)
#     make clean (effacement des fichiers objets et de l'executable)
#
 

CC       = g++-14
LD       = g++-14
CFLAGS   = -O3 
LIBS     = -L/opt/homebrew/opt/lapack/lib -llapack  ../libdvr.a

#CC       = pgCC
#LD       = pgCC
#CFLAGS    = -tp p7-64 -fastsse -O4  -pgf90libs 
#LIBS      = -L${LIBDIR} -L/home/cyril_falvo/opt/ARPACK -larpack_LINUX -lacml -lfftw3
#

OBJ_DIR = obj
SRC_DIR = src

PROG     = ex_morse.x
LISTEOBJ = \
	${OBJ_DIR}/ex_morse.o

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
	$(LD) -o $(PROG) $(CFLAGS) $(LISTEOBJ) $(LIBS)
#
# effacement des fichiers intermediaires
#
clean :
	rm -f $(PROG) $(LISTEOBJ)

