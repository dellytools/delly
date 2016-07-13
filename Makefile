DEBUG ?= 0
PARALLEL ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
SEQTK_ROOT ?= ${PWD}/src/htslib/
BOOST_ROOT ?= ${PWD}/src/modular-boost

# Flags
CXX=g++
CXXFLAGS += -isystem ${SEQTK_ROOT} -isystem ${BOOST_ROOT} -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS -fno-strict-aliasing
LDFLAGS += -L${SEQTK_ROOT} -L${BOOST_ROOT}/stage/lib -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time 

# Additional flags for release/debug
ifeq (${PARALLEL}, 1)
	CXXFLAGS += -fopenmp -DOPENMP
else
	CXXFLAGS += -DNOPENMP
endif

# Additional flags for release/debug
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz
else
	LDFLAGS += -lhts -lz -Wl,-rpath,${SEQTK_ROOT},-rpath,${BOOST_ROOT}/stage/lib
endif
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -DNDEBUG
endif


# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
BOOSTSOURCES = $(wildcard src/modular-boost/libs/iostreams/include/boost/iostreams/*.hpp)
DELLYSOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
TARGETS = .htslib .bcftools .boost src/delly src/extract src/cov src/iover src/stats src/annotate

all:   	$(TARGETS)

.htslib: $(HTSLIBSOURCES)
	cd src/htslib && make && make lib-static && cd ../../ && touch .htslib

.bcftools: $(HTSLIBSOURCES)
	cd src/bcftools && make all && cd ../../ && touch .bcftools

.boost: $(BOOSTSOURCES)
	cd src/modular-boost && ./bootstrap.sh --prefix=${PWD}/src/modular-boost --without-icu --with-libraries=iostreams,filesystem,system,program_options,date_time && ./b2 && ./b2 headers && cd ../../ && touch .boost

src/delly: .htslib .bcftools .boost $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/extract: .htslib .bcftools .boost $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/cov: .htslib .bcftools .boost $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/iover: .htslib .bcftools .boost $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/stats: .htslib .bcftools .boost $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

src/annotate: .htslib .bcftools .boost $(DELLYSOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

clean:
	cd src/htslib && make clean
	cd src/modular-boost && ./b2 --clean-all
	rm -f $(TARGETS) $(TARGETS:=.o) .htslib .boost .bcftools
