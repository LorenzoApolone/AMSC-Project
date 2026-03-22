CXX = mpic++
CXXFLAGS = -std=c++17 -O3 -Wall -pedantic -I..

TEST_CPSO = test_cpso
TEST_CPSO_PARALLEL = test_cpso_parallel

# local cpso source files
CPSO_SRCS = CPSOBase.cpp CPSOParallel.cpp CPSOSerial.cpp ContextVector.cpp SubSwarm.cpp SubSwarmTopologyFactory.cpp SwarmMetrics.cpp
CPSO_OBJS = $(CPSO_SRCS:.cpp=.o)

# external dependencies from main src/
EXTRA_SRCS = ../interfaces.cpp ../topology/create_network.cpp
EXTRA_OBJS = $(EXTRA_SRCS:.cpp=.o)

.PHONY: all clean

all: $(TEST_CPSO) $(TEST_CPSO_PARALLEL)

# Rule for building the serial test
$(TEST_CPSO): ../test_cpso.o $(EXTRA_OBJS) $(CPSO_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Rule for building the parallel test
$(TEST_CPSO_PARALLEL): ../test_cpso_parallel.o $(EXTRA_OBJS) $(CPSO_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Generic pattern rule for object files compilation
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean up build artifacts locally and contextually
clean:
	rm -f *.o $(TEST_CPSO) $(TEST_CPSO_PARALLEL)
	rm -f ../test_cpso.o ../test_cpso_parallel.o $(EXTRA_OBJS)
