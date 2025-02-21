DEBUG ?= 0
PARALLEL ?= 0
STATIC ?= 0

OS=$(shell uname)

# Submodules
PWD = $(shell pwd)
EBROOTHTSLIB ?= ${PWD}/src/htslib/

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir ?= $(exec_prefix)/bin

# Flags
CXX ?= g++
CXXFLAGS += -std=c++17 -isystem ${EBROOTHTSLIB} -pedantic -W -Wall -Wno-unknown-pragmas -D__STDC_LIMIT_MACROS -fno-strict-aliasing -fpermissive
LDFLAGS += -L${EBROOTHTSLIB} -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time

ifeq ($(OS), Darwin)
    # Skip dependency checks if target is 'clean' or 'distclean'
    ifneq ($(MAKECMDGOALS),clean)
    ifneq ($(MAKECMDGOALS),distclean)
        $(info macOS Dependency Checks...)

        ifeq ($(STATIC), 1)
            $(error Static compilation is not currently supported on macOS.)
        endif

        $(info   - checking for brew...)
        HAS_BREW := $(shell command -v brew >/dev/null 2>&1 && echo "yes" || echo "no")
        ifeq ($(HAS_BREW), no)
            $(error Homebrew is required on macOS but not found. Please install it from https://brew.sh/)
        endif

        BREW_DEPENDENCIES := autoconf automake boost llvm libomp
    
        BREW_PREFIX := $(shell brew --prefix)

        BREW_AUTOCONF_PREFIX := ${BREW_PREFIX}/opt/autoconf
        BREW_AUTOMAKE_PREFIX := ${BREW_PREFIX}/opt/automake
        BREW_BOOST_PREFIX := ${BREW_PREFIX}/opt/boost
        BREW_LLVM_PREFIX := ${BREW_PREFIX}/opt/llvm
        BREW_LIBOMP_PREFIX := ${BREW_PREFIX}/opt/libomp
        $(info   - checking for dependencies...)
        MISSING := \
            $(if $(wildcard ${BREW_AUTOCONF_PREFIX}),,autoconf) \
            $(if $(wildcard ${BREW_AUTOMAKE_PREFIX}),,automake) \
            $(if $(wildcard ${BREW_BOOST_PREFIX}),,boost) \
            $(if $(wildcard ${BREW_LLVM_PREFIX}),,llvm) \
            $(if $(wildcard ${BREW_LIBOMP_PREFIX}),,libomp)
        MISSING :=$(strip $(MISSING))
        ifneq (${MISSING},)
            $(error Install missing dependencies `brew install ${MISSING}`)
        endif

        CC = ${BREW_LLVM_PREFIX}/bin/clang
        CXX = ${BREW_LLVM_PREFIX}/bin/clang++

        CXXFLAGS += -I${BREW_BOOST_PREFIX}/include -I${BREW_LLVM_PREFIX}/include -I${BREW_LIBOMP_PREFIX}/include
        LDFLAGS += -L${BREW_BOOST_PREFIX}/lib -L${BREW_LIBOMP_PREFIX}/lib -L${BREW_LLVM_PREFIX}/lib -L${BREW_LLVM_PREFIX}/lib/c++ -L${BREW_LLVM_PREFIX}/lib/unwind -lunwind

        SDKROOT=$(shell xcrun --sdk macosx --show-sdk-path)
        SDKROOT_PREFIX = SDKROOT=$(SDKROOT)
    endif
    endif
else
    SDKROOT_PREFIX =
endif

# Flags for parallel computation
ifeq (${PARALLEL}, 1)
	CXXFLAGS += -fopenmp -DOPENMP
else
	CXXFLAGS += -DNOPENMP
endif

# Flags for static compile
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz -llzma -lbz2 -ldeflate
else
	LDFLAGS += -lhts -lz -llzma -lbz2 -Wl,-rpath,${EBROOTHTSLIB}
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
ifeq (${EBROOTHTSLIB}, ${PWD}/src/htslib/)
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
	$(SDKROOT_PREFIX) $(CXX) $(CXXFLAGS) $@.cpp src/edlib.cpp -o $@ $(LDFLAGS)

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	if [ -r src/htslib/Makefile ]; then cd src/htslib && $(MAKE) clean; fi
	rm -f $(TARGETS) $(TARGETS:=.o) ${SUBMODULES}

distclean: clean
	rm -f ${BUILT_PROGRAMS}

.PHONY: clean distclean install all
