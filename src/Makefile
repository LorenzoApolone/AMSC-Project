# Select the correct compiler
#CXX      ?= g++
CXX       = mpic++
CXXFLAGS ?= -std=c++17
CPPFLAGS ?= -O3 -Wall -pedantic -I.
LDLIBS   ?= 
LINK.o := $(LINK.cc) # Use C++ linker.

DEPEND = make.dep

# Define executables
EXEC_SERIAL = main_serial
EXEC_PARALLEL = main_parallel
EXEC_DPSO = main_dpso

.PHONY: all serial parallel dpso clean distclean

all: $(DEPEND) serial parallel dpso

serial: $(EXEC_SERIAL)

parallel: $(EXEC_PARALLEL)

dpso: $(EXEC_DPSO) main_dpso_serial

$(EXEC_SERIAL): pso/main_serial.o pso/pso_serial.o interfaces/interfaces.o dpso/dpso.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

$(EXEC_PARALLEL): pso/main_parallel.o pso/pso_mpi.o interfaces/interfaces.o dpso/dpso.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

main_dpso: dpso/main_dpso.o dpso/dpso.o interfaces/interfaces.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

main_dpso_serial: dpso/main_dpso_serial.o dpso/dpso_serial.o interfaces/interfaces_serial.o
	g++ $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

dpso/main_dpso_serial.o: dpso/main_dpso_serial.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

dpso/dpso_serial.o: dpso/dpso_serial.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

interfaces/interfaces_serial.o: interfaces/interfaces.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

interfaces/%.o: interfaces/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

pso/%.o: pso/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

dpso/%.o: dpso/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:
	$(RM) $(DEPEND)
	$(RM) *.o interfaces/*.o pso/*.o dpso/*.o

distclean: clean
	$(RM) $(EXEC_SERIAL) $(EXEC_PARALLEL) $(EXEC_DPSO) main_dpso_serial
	$(RM) *.csv *.out *.bak *~ run_parameters.json

$(DEPEND): $(wildcard *.cpp interfaces/*.cpp pso/*.cpp dpso/*.cpp)
	@$(RM) $(DEPEND)
	@for file in $(wildcard *.cpp interfaces/*.cpp pso/*.cpp dpso/*.cpp); \
	do \
	  $(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM -MT "$${file%.cpp}.o" $${file} >> $(DEPEND); \
	done

-include $(DEPEND)
