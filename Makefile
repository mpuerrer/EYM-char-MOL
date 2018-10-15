CXX = g++
CXXFLAGS = -Wall -pedantic -O3 -fomit-frame-pointer
SRC = EYM_coupled_cmp.cpp EYM_coupled_cmp_dd.cpp
INC = -I/usr/local/include
LIBS =
LIBS_DD = -lqd
SRC = EYM_coupled_cmp.cpp EYM_coupled_cmp_dd.cpp

.SUFFIXES : .o .cpp
.cpp.o :
	$(CXX) $(CXXFLAGS) $(INC) -c $<

.PHONY: all clean depend

all: EYM_coupled_cmp EYM_coupled_cmp_dd

EYM_coupled_cmp: EYM_coupled_cmp.o
	$(CXX) -o EYM_coupled_cmp EYM_coupled_cmp.o $(LIBS)

EYM_coupled_cmp_dd: EYM_coupled_cmp_dd.o
	$(CXX) -o EYM_coupled_cmp_dd EYM_coupled_cmp_dd.o $(LIBS_DD)


clean:
	rm -f *.o *.dat EYM_coupled_cmp EYM_coupled_cmp_dd

depend:
	makedepend -Y $(SRC)
# DO NOT DELETE
