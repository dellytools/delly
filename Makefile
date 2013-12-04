DEBUG ?= 0

# External libraries
BOOST=/g/solexa/bin/software/boost_1_53_0
BAMTOOLS=/g/solexa/bin/software/bamtools-2.3.0/
KSEQ=/g/solexa/bin/software/kseq/

# Flags
CXX=g++
CXXFLAGS += -isystem ${BOOST}/include -isystem ${BAMTOOLS}/include -isystem ${KSEQ}
LDFLAGS += -L${BOOST}/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time -L${BAMTOOLS}/lib -lbamtools -lz

# Additional flags for release/debug
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -pedantic -W -Wall -DDEBUG
	LDFLAGS += -lprofiler
else
	CXXFLAGS += -O9 -pedantic -W -Wall -DNDEBUG
	LDFLAGS += --static
endif

### Memory debugging: Valgrind
## Set the library path: export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${BOOST}/lib:${BAMTOOLS}/lib
## Run valgrind (1): valgrind --leak-check=yes --track-origins=yes --num-callers=30 --error-limit=no --track-fds=yes ...
## Run valgrind (2): valgrind --tool=massif ...
### Code checker: Cppcheck
## cppcheck -I${BAMTOOLS}/include -I${KSEQ} -I${BOOST}/include/ --enable=style,performance,portability,information,unusedFunction delly/delly.cpp
### Profiling: gperftools 
## CPUPROFILE=delly.prof ./delly ...
## google-pprof --gv src/delly delly.prof

TARGETS = src/delly src/iover src/cov src/spancov src/iMerge src/extract 

all:   	$(TARGETS)

src/delly:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/iover:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/iMerge:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/extract:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/cov:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/spancov:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGETS) $(TARGETS:=.o)
