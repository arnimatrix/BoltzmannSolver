VERSION=1.0

CXX 		 = g++-7
CXXDEBUG = $(CXX) -g -D DEBUG
CXXFLAGS = -Wall -Wextra -Wpedantic -std=c++17 -O3 -fPIC

SCRIPTDIR = script
BINDIR    = bin

# Paths to headers (-I...)
IPATH = 

# Paths to libraries (-L...)
LPATH =

# Libraries
LIBS  = -lmarty
 
SCRIPTS  = $(wildcard $(SCRIPTDIR)/*.cpp)
BIN      = $(SCRIPTS:$(SCRIPTDIR)/%.cpp=$(BINDIR)/%.x) 
BINDEBUG = $(SCRIPTS:$(SCRIPTDIR)/%.cpp=$(BINDIR)/%_debug.x) 
	
release: $(BIN)
debug: $(BINDEBUG)

all: $(BIN) $(BINDEBUG)

# Executables
$(BINDIR)/%.x: $(SCRIPTDIR)/%.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(IPATH) $(LPATH) $(LIBS)

# Debugging version
$(OBJDIR)/%_debug.x: $(SCRIPTDIR)/%.cpp
	$(CXXDEBUG) $(CXXFLAGS) $< -o $@ $(IPATH) $(LPATH) $(LIBS)

.PHONY: clean
clean:
	-rm -f $(BINDIR)/*.x;
