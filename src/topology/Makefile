CXX = mpic++
CXXFLAGS = -O2 -std=gnu++17

TARGET = topology_parallel
TARGET_SERIAL = topology_serial
SRCS = main_topology.cpp create_network.cpp confront.cpp pso_serial_topology.cpp pso_topology.cpp ../interfaces/interfaces.cpp ../pso/pso_mpi.cpp ../pso/pso_serial.cpp
SRCS_SERIAL = main_topology_serial.cpp create_network.cpp confront.cpp pso_serial_topology.cpp pso_topology.cpp ../interfaces/interfaces.cpp ../pso/pso_mpi.cpp ../pso/pso_serial.cpp

all: $(TARGET) $(TARGET_SERIAL)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

$(TARGET_SERIAL): $(SRCS_SERIAL)
	$(CXX) $(CXXFLAGS) $(SRCS_SERIAL) -o $(TARGET_SERIAL)

clean:
	rm -f $(TARGET) $(TARGET_SERIAL)
