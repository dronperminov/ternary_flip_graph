CXX = g++
FLAGS = -Wall -O3 -std=c++17 -fopenmp
OBJECTS = src/entities/arg_parser.o src/entities/flip_set.o src/utils.o src/schemes/base_scheme.o

all: ternary_flip_graph ternary_meta_flip_graph find_alternative_schemes

ternary_flip_graph: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) main_flip_graph.cpp -o ternary_flip_graph

ternary_meta_flip_graph: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) main_meta_flip_graph.cpp -o ternary_meta_flip_graph

find_alternative_schemes: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) find_alternative_schemes.cpp -o find_alternative_schemes

%.o: %.cpp
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm -rf $(OBJECTS) ternary_flip_graph ternary_meta_flip_graph find_alternative_schemes
