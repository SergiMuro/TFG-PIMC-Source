# to compile with icc run $> module load toolchains/intel_mkl_ompi
CCOMPILER = gcc 
#CCOMPILER = mpicc
FCOMPILER = cc
#FLAGS = -g
#FLAGS = -O2 -static-intel
#FLAGS = -O3 -openmp
# use -g for valgrind
#FLAGS = -O3 -qopenmp
#FLAGS = -O3 -fopenmp
#FLAGS = -O3
#FLAGS = -O0 -g
# -Wall -ffast-math
#FLAGS = -g
# for intel debugger use idb and compile with flag:
#FLAGS = -debug -O0
# Xean (Skylake) -O3 -qopenmp -xCORE-AVX512 . 
# Xeon Phi (KNL) -xMIC-AVX512
#
MPI_LIBS=/usr/lib/lam/
MPI_INCLUDE_DIR=/usr/include/lam/

LIBRARIES = -L$(MPI_LIBS) -lm
#LIBRARIES = -L$(MPI_LIBS) -lmpich -lm
#ifdef __linux__
  GRX20 = c:/djgpp/lib/libgrx20.a
  BCC = c:/djgpp/lib/libbcc.a
  BCC2 = c:/djgpp/lib/dj/libbcc2.a
#endif
  GRX20 = /usr/lib/libgrx20.a
  BCC = libbcc.a
  BCC2 = /usr/lib/libbcc2.a
#else

SRC = $(*.c)
OBJ = dmc.o display.o gencoord.o memory.o vmc.o randnorm.o \
      rw.o trial.o utils.o quantities.o parallel.o crystal.o \
      spline.o optimiz.o mymath.o ewald.o main.o pimc.o tvmc.o pigs.o libigor.o

all: MC

#ifdef __linux__
  OUT = MC
#else
  OUT = MC.exe
#endif

MC : $(OBJ)
	$(CCOMPILER) $(OBJ) $(FLAGS) $(LIBRARIES) -o MC 

.o:
	$(CCOMPILER) -I$(MPI_INCLUDE_DIR) $(FLAGS) $*.c -o $*.o -Wall -c

$(OBJ): main.c
	$(CCOMPILER) -I$(MPI_INCLUDE_DIR) $(FLAGS) $*.c -o $*.o -Wall -c

#ran2.o: ran2.f
#	$(FCOMPILER) ran2.f -c
        
clean :
	rm *.o	

