
# Makefile

CC = icpc

MKLROOT = /opt/intel/composer_xe_2013_sp1/mkl

CFLAGS = -g -debug -O2 -openmp -I$(MKLROOT)/include -mkl

BOOST_INC = -I/opt/local/include/

OBJ  =   obj/main.o  obj/matrix.o  

LIBS = -lm

HEADS =  include/matrix.hpp
BIN  =   DMRG

RM = rm  -f

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
