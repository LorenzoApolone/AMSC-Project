CXX = mpic++
CXXFLAGS = -O2 -std=gnu++17

TARGET = a
SRCS = main_topology.cpp create_network.cpp confront.cpp ../interfaces.cpp ../functions.cpp ../pso_mpi.cpp

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET)

