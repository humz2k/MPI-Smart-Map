BUILD ?= build
SOURCE ?= src
INCLUDE ?= include

MPI_CXX ?= mpicxx

INCLUDES ?= -I$(INCLUDE)

SOURCES := $(shell find $(SOURCE) -name '*.cpp')
OBJECTS := $(SOURCES:%.cpp=%.o)
OUTPUTS := $(OBJECTS:%=$(BUILD)/%)

HEADERS := $(shell find $(INCLUDE) -name '*.hpp')

driver: $(OUTPUTS)
	$(MPI_CXX) $^ -o $@

$(BUILD)/%.o: %.cpp $(HEADERS)
	mkdir -p $(@D)
	$(MPI_CXX) -c $< $(INCLUDES) -o $@

.PHONY: test
test: driver
	mpirun -n 8 driver

.PHONY: clean
clean:
	rm -rf build
	rm driver