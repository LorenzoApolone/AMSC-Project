CXX = mpic++
CXXFLAGS = -O2 -std=gnu++17

TARGET = a

SRCS = main_small.cpp create_network.cpp ../interfaces.cpp ../functions.cpp

all: $(TARGET)

$(TARGET):
	$(CXX) $(CXXFLAGS) $(SRCS) -o $(TARGET)

clean:
	rm -f $(TARGET)

