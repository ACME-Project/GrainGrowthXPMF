# Makefile

# includes
incdir = $(MMSP_PATH)/include
utildir = $(MMSP_PATH)/utility
algodir = $(MMSP_PATH)/algorithms

# compilers/flags
compiler = g++ -O3 -Wall
pcompiler = mpic++ -O3 -Wall
flags = -I$(incdir) -I$(algodir) -I$(utildir)

# IBM compiler for AMOS
BG_XL = /bgsys/drivers/ppcfloor/comm/xl
BG_INC = -I$(BG_XL)/include
BG_LIB = -L$(BG_XL)/lib
qcompiler = $(BG_XL)/bin/mpixlcxx_r -O3 -qflag=w -qstrict -qmaxmem=-1
qflags = $(BG_INC) $(BG_LIB) $(flags) -I/bgsys/apps/CCNI/zlib/zlib-1.2.7/include -L/bgsys/apps/CCNI/zlib/zlib-1.2.7/lib

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp \
       $(incdir)/MMSP.sparse.hpp \

# the program
sparse: sp-xmpf.cpp $(core)
	$(pcompiler) $(flags) -include mpi.h $< -o sparse.out -lz

bgq: sp-xmpf.cpp $(core)
	$(qcompiler) -DBGQ $(qflags) $< -o q_GG.out -lz

clean:
	rm -rf q_GG.out sparse.out
