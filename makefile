CXX = g++
FLAGS = -Wall -O3 -std=c++17 -fopenmp
OBJECTS = src/entities/arg_parser.o src/entities/flip_set.o

all: ternary_flip_graph

ternary_flip_graph: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) main.cpp -o ternary_flip_graph

%.o: %.cpp
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm -rf src/*/*.o ternary_flip_graph
