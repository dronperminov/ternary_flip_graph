CXX = g++
FLAGS = -Wall -O3 -std=c++17 -fopenmp
ALGEBRA_OBJECTS = src/algebra/fraction.o src/algebra/matrix.o src/algebra/binary_matrix.o src/algebra/binary_solver.o src/algebra/mod3_solver.o
ENTITIES_OBJECTS = src/entities/arg_parser.o src/entities/flip_set.o src/entities/ranks.o src/entities/flip_parameters.o src/entities/meta_parameters.o
LIFT_OBJECTS = src/lift/binary_lifter.o src/lift/mod3_lifter.o
SCHEMES_OBJECTS = src/schemes/base_scheme.o src/schemes/fractional_scheme.o
OBJECTS = $(ALGEBRA_OBJECTS) ${ENTITIES_OBJECTS} $(LIFT_OBJECTS) $(SCHEMES_OBJECTS) src/utils.o

all: flip_graph meta_flip_graph complexity_minimizer find_alternative_schemes lift

flip_graph: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) flip_graph.cpp -o flip_graph

meta_flip_graph: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) meta_flip_graph.cpp -o meta_flip_graph

complexity_minimizer: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) complexity_minimizer.cpp -o complexity_minimizer

find_alternative_schemes: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) find_alternative_schemes.cpp -o find_alternative_schemes

validate_schemes: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) validate_schemes.cpp -o validate_schemes

lift: $(OBJECTS)
	$(CXX) $(FLAGS) $(OBJECTS) lift.cpp -o lift

%.o: %.cpp
	$(CXX) $(FLAGS) -c $< -o $@

clean:
	rm -rf $(OBJECTS) flip_graph meta_flip_graph complexity_minimizer find_alternative_schemes
