CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

TARGET = solver
SRCS = main.cpp Params.cpp Solver.cpp
OBJS = $(SRCS:.cpp=.o)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

main.o: main.cpp Params.h Field.h Index.h Solver.h
	$(CXX) $(CXXFLAGS) -c main.cpp

Params.o: Params.cpp Params.h
	$(CXX) $(CXXFLAGS) -c Params.cpp

Solver.o: Solver.cpp Solver.h Params.h Field.h Index.h
	$(CXX) $(CXXFLAGS) -c Solver.cpp

clean:
	rm -f $(OBJS) $(TARGET)