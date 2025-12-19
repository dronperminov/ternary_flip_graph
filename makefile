CXX = g++
FLAGS = -Wall -O3 -std=c++17 -fopenmp
OBJECTS = src/entities/arg_parser.o src/entities/flip_set.o src/utils.o

all: ternary_flip_graph ternary_meta_flip_graph

ternary_flip_graph: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) main_flip_graph.cpp -o ternary_flip_graph

ternary_meta_flip_graph: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) main_meta_flip_graph.cpp -o ternary_meta_flip_graph

%.o: %.cpp
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	del /F /S /Q *.o ternary_flip_graph ternary_meta_flip_graph
