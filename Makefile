# Density Matrix Renormalization Group (DMRG) Code
# Authors J. Szekely, B. Ashwell and S. Parker
# Makefile

CC = g++

#MKLROOT = /opt/intel/composer_xe_2013_sp1/mkl

CFLAGS = -O2 -I$(MKLROOT)/include -Wall -Wno-sign-compare -Wno-unused-function -Werror -std=c++11 -openmp

BOOST_INC = -I/opt/local/include/

OBJ  =   obj/main.o  obj/matrix.o obj/dmrg.o obj/davidson.o

LIBS =  -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_sequential -lpthread -lm

HEADS =  include/matrix.hpp include/utilities.hpp include/dmrg.hpp include/davidson.hpp
BIN  =   DMRG

RM = rm -f

.PHONY: all all-before all-after clean clean-custom
all: all-before $(BIN) all-after

clean: clean-custom
	${RM} $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CC) $(CFLAGS) -o $(BIN) -I./include $(OBJ) $(BOOST_INC) $(LIBS)

obj/main.o: src/main.cpp
	$(CC) $(CFLAGS) -c  src/main.cpp -o obj/main.o $(BOOST_INC) -I./include

obj/matrix.o: src/matrix.cpp
	$(CC) $(CFLAGS)  -c src/matrix.cpp  -o obj/matrix.o  -I./include

obj/dmrg.o: src/dmrg.cpp
	$(CC) $(CFLAGS)  -c src/dmrg.cpp  -o obj/dmrg.o  -I./include

obj/davidson.o: src/davidson.cpp
	$(CC) $(CFLAGS)  -c src/davidson.cpp  -o obj/davidson.o  -I./include
