
CHPL = chpl
FLAGS = --fast

PAR_FLAGS =  
MUL_FLAGS =

HEAT_PAR = heat
HEAT_MUL = heat_multi

heat_SOURCES = main.chpl common.chpl compute.chpl

EXAMPLE_PAR = example_chpl
EXAMPLE_MUL = example_chpl_multi

example_SOURCES = example.chpl

HELLO_LOC = hello_loc
HELLO_DIST = hello_dist

hello_SOURCES = hello.chpl

all: $(CONF) $(heat_SOURCES) $(example_sources) $(hello_SOURCES) $(HEAT_SEQ) $(HEAT_PAR) $(HEAT_MUL) $(EXAMPLE_SEQ) $(EXAMPLE_PAR) $(EXAMPLE_MUL) $(HELLO_LOC) $(HELLO_DIST)

$(HEAT_PAR): $(heat_SOURCES)
	$(CHPL) -o $(HEAT_PAR) $(FLAGS) $(PAR_FLAGS) $(heat_SOURCES)

$(HEAT_MUL): $(heat_SOURCES)
	CHPL_COMM_SUBSTRATE=ibv CHPL_COMM=gasnet $(CHPL) -o $(HEAT_MUL) $(FLAGS) $(MUL_FLAGS) $(heat_SOURCES)

$(EXAMPLE_PAR): $(example_SOURCES)
	$(CHPL) -o $(EXAMPLE_PAR) $(FLAGS) $(PAR_FLAGS) $(example_SOURCES)

$(EXAMPLE_MUL): $(example_SOURCES)
	CHPL_COMM=gasnet $(CHPL) -o $(EXAMPLE_MUL) $(FLAGS) $(MUL_FLAGS) $(example_SOURCES)

$(HELLO_LOC): $(hello_SOURCES)
	$(CHPL) -o $(HELLO_LOC) $(FLAGS) $(hello_SOURCES)

$(HELLO_DIST): $(hello_SOURCES)
	CHPL_COMM=gasnet $(CHPL) -o $(HELLO_DIST) $(FLAGS) $(hello_SOURCES)

clean: 
	rm -f $(EXAMPLE_PAR) $(EXAMPLE_MUL) $(HEAT_PAR) $(HEAT_MUL) $(HEAT_MUL)_real $(EXAMPLE_MUL)_real $(HELLO_LOC){,_real} $(HELLO_DIST){,_real}


