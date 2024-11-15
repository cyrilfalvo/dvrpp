#
#  Makefile 
#     make all   (construction de l'executable)
#     make clean (effacement des fichiers objets et de l'executable)
#

CC       = g++-14
AR       = ar
CFLAGS   = -O2 

##CFLAGS   = -O3 -D_ARPACK_

OBJ_DIR = obj
SRC_DIR = src

LIBDVR = libdvr.a
LISTEOBJ = \
        ${OBJ_DIR}/BasisSet.o \
        ${OBJ_DIR}/BasisHerm.o \
        ${OBJ_DIR}/BasisSine.o \
        ${OBJ_DIR}/BasisExp.o \
        ${OBJ_DIR}/DVRpp.o \
	${OBJ_DIR}/DVRpp_wrapper.o\

#
# compilation
#
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

#
# dependances
#
all       : $(LIBDVR)
#
# edition de lien
#
$(LIBDVR) : $(LISTEOBJ)
	$(AR) rcs $(LIBDVR) $(LISTEOBJ) 


# effacement des fichiers intermediaires
#
clean :
	rm -f $(LIBDVR) $(LISTEOBJ)

