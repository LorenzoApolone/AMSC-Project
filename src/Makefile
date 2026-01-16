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

# Common source files (excluding main files and algorithm implementations)
COMMON_SRCS := $(filter-out main_serial.cpp main_parallel.cpp pso_serial.cpp pso_mpi.cpp, $(wildcard *.cpp))
COMMON_OBJS = $(COMMON_SRCS:.cpp=.o)

# Source files for serial version
SERIAL_SRCS = $(COMMON_SRCS) main_serial.cpp pso_serial.cpp
SERIAL_OBJS = $(SERIAL_SRCS:.cpp=.o)

# Source files for parallel version  
PARALLEL_SRCS = $(COMMON_SRCS) main_parallel.cpp pso_mpi.cpp
PARALLEL_OBJS = $(PARALLEL_SRCS:.cpp=.o)

.PHONY = all serial parallel clean distclean $(DEPEND)

all: $(DEPEND) serial parallel

serial: $(EXEC_SERIAL)

parallel: $(EXEC_PARALLEL)

$(EXEC_SERIAL): $(SERIAL_OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

$(EXEC_PARALLEL): $(PARALLEL_OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

clean:
	$(RM) $(DEPEND)
	$(RM) *.o

distclean: clean
	$(RM) $(EXEC_SERIAL) $(EXEC_PARALLEL)
	$(RM) *.csv *.out *.bak *~ run_parameters.json

$(DEPEND): $(wildcard *.cpp)
	@$(RM) $(DEPEND)
	@for file in $(wildcard *.cpp); \
	do \
	  $(CXX) $(CPPFLAGS) $(CXXFLAGS) -MM $${file} >> $(DEPEND); \
	done

-include $(DEPEND)
