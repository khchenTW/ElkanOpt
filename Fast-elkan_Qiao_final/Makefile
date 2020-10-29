# compile up to 8 things at once
MAKEFLAGS = -j 8

# Compiler options
CPPFLAGS = -Werror -pedantic
CPPFLAGS += -fPIC
CPPFLAGS += -g
CPPFLAGS += -O3
CPPFLAGS += -Wno-long-long
CPPFLAGS += -std=c++11

# Verify algorithm correctness while debugging
#CPPFLAGS += -DDEBUG
#CPPFLAGS += -DVERIFY_ASSIGNMENTS

# To use pthreads, uncomment both lines below.
#CPPFLAGS += -DUSE_THREADS
#LDFLAGS += -lpthread

# Monitor internal algorithm effectiveness
#CPPFLAGS += -DCOUNT_DISTANCES
#CPPFLAGS += -DMONITOR_ACCURACY

# Enable code profiling
#CPPFLAGS += -pg

# Source files for libkmeans, excluding files for Python module

LIBSRC = $(shell ls ./src/*.cpp | grep -v driver | grep -v fastkmeans)
LIBOBJS = $(LIBSRC:.cpp =.o)

KMEANSLIBRARY = ./bin/libkmeans.a

all: driver-standalone

$(KMEANSLIBRARY): $(LIBOBJS)
	ar -cr $(KMEANSLIBRARY) $(LIBOBJS)

driver-standalone: $(KMEANSLIBRARY) driver-standalone.o
	g++ -L . $(CPPFLAGS) $(LDFLAGS) ./bin/driver-standalone.o -o driver-standalone -lkmeans


.PHONY: clean all

clean:
	rm -f $(KMEANSLIBRARY) driver-standalone *.o gmon.out

