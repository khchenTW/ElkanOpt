# compile up to 8 things at once
MAKEFLAGS = -j 8

# Compiler and preprocessor flags
CC:= g++
CPPFLAGS = -Werror -pedantic -Iinclude/
CPPFLAGS += -fPIC
CPPFLAGS += -g
CPPFLAGS += -O3
CPPFLAGS += -Wno-long-long
CPPFLAGS += -std=c++11
LDLIBS:= -llibkmeans
LDFLAGS:= -Llib

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

SRC_DIR:=src
OBJ_DIR:=obj
BIN_DIR:=bin

EXE:=$(BIN_DIR)/driver-standalone
SRC:=$(wildcard $(SRC_DIR)/*.cpp | grep -v driver | grep -v fastkmeans)
OBJ:=$(SRC:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)
# Source files for libkmeans, excluding files for Python module


KMEANSLIBRARY = ./libkmeans.a

.PHONY: clean all

all: $(EXE)

$(KMEANSLIBRARY): $(OBJ)
	ar -cr $(KMEANSLIBRARY) $(OBJ)

$(EXE): $(OBJ) | $(BIN_DIR) $(KMEANSLIBRARY)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CC) $(CPPFLAGS) -c -o $@ $<

$(BIN_DIR):
	mkdir -p $@

$(OBJ_DIR):
	mkdir -p $@

python-module: $(KMEANSLIBRARY)
	cd python-bindings && python3 setup.py build_ext --inplace

clean:
	rm -rv $(KMEANSLIBRARY) $(OBJ_DIR) $(BIN_DIR)

-include $(OBJ:.o=.d)
