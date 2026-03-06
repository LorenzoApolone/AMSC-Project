# Select the correct compiler
#CXX      ?= g++
CXX       = mpic++
CXXFLAGS ?= -std=c++17
CPPFLAGS ?= -O3 -Wall -pedantic -I.
LDLIBS   ?= 
LINK.o := $(LINK.cc) # Use C++ linker.

DEPEND = make.dep

# Define both executables
EXEC_SERIAL = main_serial
EXEC_PARALLEL = main_parallel
EXEC_DPSO = main_dpso

# Common source files (excluding main files and algorithm implementations)
COMMON_SRCS := $(filter-out main_serial.cpp main_parallel.cpp pso_serial.cpp pso_mpi.cpp, $(wildcard *.cpp))
COMMON_OBJS = $(COMMON_SRCS:.cpp=.o)

# Source files for serial version
SERIAL_SRCS = $(COMMON_SRCS) main_serial.cpp pso_serial.cpp
SERIAL_OBJS = $(SERIAL_SRCS:.cpp=.o)

# Source files for parallel version  
PARALLEL_SRCS = $(COMMON_SRCS) main_parallel.cpp pso_mpi.cpp
PARALLEL_OBJS = $(PARALLEL_SRCS:.cpp=.o)

# Aggiungi i file della cartella dpso
DPSO_SRCS := dpso/dpso.cpp
DPSO_OBJS := $(DPSO_SRCS:.cpp=.o)

# Aggiorna le variabili dei sorgenti
COMMON_SRCS := $(filter-out main_serial.cpp main_parallel.cpp pso_serial.cpp pso_mpi.cpp, $(wildcard *.cpp)) $(DPSO_SRCS)
COMMON_OBJS := $(COMMON_SRCS:.cpp=.o)

# Aggiorna le variabili degli oggetti per seriale e parallelo
SERIAL_SRCS = $(COMMON_SRCS) main_serial.cpp pso_serial.cpp
SERIAL_OBJS = $(SERIAL_SRCS:.cpp=.o)

PARALLEL_SRCS = $(COMMON_SRCS) main_parallel.cpp pso_mpi.cpp
PARALLEL_OBJS = $(PARALLEL_SRCS:.cpp=.o)

.PHONY = all serial parallel clean distclean $(DEPEND)

all: $(DEPEND) serial parallel

serial: $(EXEC_SERIAL)

parallel: $(EXEC_PARALLEL)

$(EXEC_SERIAL): main_serial.o pso_serial.o functions.o interfaces.o dpso/dpso.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

$(EXEC_PARALLEL): main_parallel.o pso_mpi.o functions.o interfaces.o dpso/dpso.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

main_dpso: dpso/main_dpso.o dpso/dpso.o functions.o interfaces.o
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

main_dpso_serial: dpso/main_dpso_serial.o dpso/dpso_serial.o functions_serial.o interfaces_serial.o
	g++ $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

dpso/main_dpso_serial.o: dpso/main_dpso_serial.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

dpso/dpso_serial.o: dpso/dpso_serial.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

interfaces_serial.o: interfaces.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

functions_serial.o: functions.cpp
	g++ $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

dpso/%.o: dpso/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:
	$(RM) $(DEPEND)
	$(RM) *.o dpso/*.o

distclean: clean
	$(RM) $(EXEC_SERIAL) $(EXEC_PARALLEL) $(EXEC_DPSO) main_dpso_serial
	$(RM) *.csv *.out *.bak *~ run_parameters.json

$(DEPEND): $(wildcard *.cpp)
	@$(RM) $(DEPEND)
	@for file in $(wildcard *.cpp); \
	do \
	  $(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM $${file} >> $(DEPEND); \
	done

-include $(DEPEND)
