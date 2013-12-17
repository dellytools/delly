DEBUG ?= 0
PARALLEL ?= 0

# External libraries
BOOST=/g/solexa/bin/software/boost_1_53_0
BAMTOOLS=/g/solexa/bin/software/bamtools-2.3.0/
KSEQ=/g/solexa/bin/software/kseq/

# Flags
CXX=g++
CXXFLAGS += -isystem ${BOOST}/include -isystem ${BAMTOOLS}/include -isystem ${KSEQ} -pedantic -W -Wall -Wno-unknown-pragmas
LDFLAGS += -L${BOOST}/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time -L${BAMTOOLS}/lib -lbamtools -lz

# Additional flags for release/debug
ifeq (${PARALLEL}, 1)
	CXXFLAGS += -fopenmp -DOPENMP
else
	CXXFLAGS += -DNOPENMP
endif

# Additional flags for release/debug
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O9 -DNDEBUG
	#LDFLAGS += --static
endif

# Targets
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
