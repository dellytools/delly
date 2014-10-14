DEBUG ?= 0
PARALLEL ?= 0

# External libraries

BOOST_ROOT ?= /g/solexa/bin/software/boost_1_53_0
BAMTOOLS_ROOT = src/bamtools
SEQTK_ROOT    = src/htslib/htslib

# Flags
CXX=g++
CXXFLAGS += -isystem ${BOOST_ROOT}/include -isystem ${BAMTOOLS_ROOT}/include -isystem ${SEQTK_ROOT} -pedantic -W -Wall -Wno-unknown-pragmas
LDFLAGS += -L${BOOST_ROOT}/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time -L${BAMTOOLS_ROOT}/lib -lbamtools -lz

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
TARGETS = htslib bamtools libbamtools.a src/delly src/extract src/cov src/iover src/stats

all:   	$(TARGETS)

htslib:
	cd src/htslib && make
bamtools:
	cd src/bamtools && mkdir -p build && cd build && cmake .. && make

libbamtools.a: bamtools
	cp src/bamtools/lib/libbamtools.a .

src/delly:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/extract:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/cov:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/iover:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/stats:
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	rm -f $(TARGETS) $(TARGETS:=.o)
