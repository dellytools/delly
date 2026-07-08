# To build against HTSlib system libraries
#
#   make HTSLIBINCDIR=/usr/include HTSLIBLIBDIR=/usr/lib all
#
PWD = $(shell pwd)
HTSLIBINCDIR ?= ${PWD}/src/htslib/
HTSLIBLIBDIR ?= ${PWD}/src/htslib/
BOOSTINCDIR ?=
BOOSTLIBDIR ?=

# Debug or static build
DEBUG ?= 0
STATIC ?= 0

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir ?= $(exec_prefix)/bin

# Flags
CXX ?= g++
CXXFLAGS += -std=c++17 -isystem ${HTSLIBINCDIR} -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS -fno-strict-aliasing -fpermissive -pthread
LDFLAGS += -L${HTSLIBLIBDIR} -lboost_iostreams -lboost_filesystem -lboost_program_options -lboost_date_time -pthread

# Boost location
ifneq (${BOOSTINCDIR},)
	CXXFLAGS += -isystem ${BOOSTINCDIR}
endif
ifneq (${BOOSTLIBDIR},)
	LDFLAGS += -L${BOOSTLIBDIR} -Wl,-rpath,${BOOSTLIBDIR}
endif

# Flags for static compile
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -lhts -lz -llzma -lbz2 -ldeflate
else
	LDFLAGS += -lhts -lz -llzma -lbz2 -Wl,-rpath,${HTSLIBLIBDIR}
endif

# Flags for debugging, profiling and releases
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG
endif
ifeq (${HTSLIBINCDIR}, ${PWD}/src/htslib/)
	SUBMODULES += .htslib
endif

# External sources
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
SOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
BUILT_PROGRAMS = src/delly
TARGETS = ${SUBMODULES} ${BUILT_PROGRAMS}

all:   	$(TARGETS)

.htslib: $(HTSLIBSOURCES)
	if [ -r src/htslib/Makefile ]; then cd src/htslib && autoreconf -i && ./configure --disable-s3 --disable-gcs --disable-libcurl --disable-plugins && $(MAKE) && $(MAKE) lib-static && cd ../../ && touch .htslib; fi

src/delly: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp src/edlib.cpp -o $@ $(LDFLAGS)

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	if [ -r src/htslib/Makefile ]; then cd src/htslib && $(MAKE) clean; fi
	rm -f $(TARGETS) $(TARGETS:=.o) ${SUBMODULES}

distclean: clean
	rm -f ${BUILT_PROGRAMS}

.PHONY: clean distclean install all
