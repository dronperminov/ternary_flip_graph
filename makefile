CXX = g++
FLAGS = -Wall -O3 -std=c++14 -fopenmp
OBJECTS = src/entities/arg_parser.o src/entities/flip_set.o src/schemes/ternary_scheme.o src/entities/flip_graph.o

all: ternary_flip_graph

ternary_flip_graph: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) main.cpp -o ternary_flip_graph

%.o: %.cpp
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm -rf src/*/*.o ternary_flip_graph
