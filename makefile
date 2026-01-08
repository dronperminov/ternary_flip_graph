CXX = g++
FLAGS = -Wall -O3 -std=c++17 -fopenmp
OBJECTS = src/algebra/fraction.o src/algebra/matrix.o src/algebra/binary_matrix.o src/entities/arg_parser.o src/entities/flip_set.o src/utils.o src/schemes/base_scheme.o

all: flip_graph meta_flip_graph complexity_minimizer find_alternative_schemes lift

flip_graph: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) flip_graph.cpp -o flip_graph

meta_flip_graph: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) meta_flip_graph.cpp -o meta_flip_graph

complexity_minimizer: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) complexity_minimizer.cpp -o complexity_minimizer

find_alternative_schemes: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) find_alternative_schemes.cpp -o find_alternative_schemes

lift: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) lift.cpp -o lift

%.o: %.cpp
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm -rf $(OBJECTS) flip_graph meta_flip_graph complexity_minimizer find_alternative_schemes
