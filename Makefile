DEBUG ?= 0
PARALLEL ?= 0

# Boost library
BOOST_ROOT ?= /g/solexa/bin/software/boost_1_53_0

# Submodules
PWD = $(shell pwd)
BAMTOOLS_ROOT = ${PWD}/src/bamtools
SEQTK_ROOT    = ${PWD}/src/htslib/htslib

# Flags
CXX=g++
CXXFLAGS += -isystem ${BOOST_ROOT}/include -isystem ${BAMTOOLS_ROOT}/include -isystem ${SEQTK_ROOT} -pedantic -W -Wall -Wno-unknown-pragmas
LDFLAGS += -L${BOOST_ROOT}/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time -L${BAMTOOLS_ROOT}/lib -lbamtools -lz -Wl,-rpath,${BAMTOOLS_ROOT}/lib,-rpath,${BOOST_ROOT}/lib

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

# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
BAMTOOLSSOURCES = $(wildcard src/bamtools/src/api/*.h) $(wildcard src/bamtools/src/api/*.cpp)
DELLYSOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
TARGETS = .htslib .bamtools src/delly src/extract src/cov src/iover src/stats

all:   	$(TARGETS)

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && cd ../../ && touch .htslib

.bamtools: $(BAMTOOLSSOURCES)
	cd src/bamtools && mkdir -p build && cd build && cmake .. && make && cd ../../../ && touch .bamtools

src/delly: .htslib .bamtools $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/extract: .htslib .bamtools $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/cov: .htslib .bamtools $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/iover: .htslib .bamtools $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/stats: .htslib .bamtools $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	cd src/bamtools/build && make clean
	cd src/htslib && make clean
	rm -f $(TARGETS) $(TARGETS:=.o)
