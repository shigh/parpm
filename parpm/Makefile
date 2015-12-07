
CXX         = mpic++
CXXFLAGS    = -g -fPIC -fopenmp
CXXINCLUDES = -I$(HOME)/local/include
CXXDEFS     = -DHAVE_CONFIG_H

LDLIBS  = -L$(HOME)/local/lib  -lfftw3_mpi -lfftw3 -lm

DEPS = fftwmpi.hpp core.hpp grid.hpp interp.hpp particle.hpp types.hpp
OBJS = fftwmpi.o core.o interp.o particle.o

pic-fft: pic-fft.o $(OBJS)
	$(CXX)  $(CXXFLAGS) -o $@ $^ $(LDLIBS)
pic-fft.o : pic-fft.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o $@ -c $<

test: test-solvers.o problem.o $(OBJS)
	$(CXX)  $(CXXFLAGS) -o $@ $^ $(LDLIBS)
problem.o: problem.hpp
run-test: test
	mpirun -n 2 ./test

test-solvers.o : test-solvers.cpp $(DEPS)
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o $@ -c $<

%.o : %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(CXXINCLUDES) -o $@ -c $<

clean:
	rm -vf *.o *.so *~
	rm -vrf build
	rm -vf test pic-fft
