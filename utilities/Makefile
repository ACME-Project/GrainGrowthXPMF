# Makefile

# includes
incdir = $(MMSP_PATH)/include
utildir = $(MMSP_PATH)/utility
algodir = $(MMSP_PATH)/algorithms

# compilers/flags
compiler = g++ -O3 -Wall
pcompiler = mpic++ -O3 -Wall
flags = -I$(incdir) -I$(algodir) -I$(utildir)

# dependencies
core = $(incdir)/MMSP.main.hpp \
       $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp \
       $(incdir)/MMSP.sparse.hpp \
       $(utildir)/tessellate.hpp

tool: tool.cpp $(core) /usr/include/IL/devil_cpp_wrapper.hpp
	$(compiler) $(flags) -I /usr/include/IL -include il.h $< -o $@ -lz -lIL -lILU -lILUT
