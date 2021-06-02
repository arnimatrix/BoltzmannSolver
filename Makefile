VERSION=1.0

CXX 		 = g++-7
CXXDEBUG = $(CXX) -g -D DEBUG
CXXFLAGS = -Wall -Wextra -Wpedantic -std=c++17 -O3 -fPIC

SCRIPTOBJDIR = script/obj
SCRIPTDIR = script
SRCDIR    = src
INCDIR    = include
OBJDIR    = obj
BINDIR    = bin

# Paths to headers (-I...)
IPATH = -I$(INCDIR)

# Paths to libraries (-L...)
LPATH =

# Libraries
LIBS  = -lmarty
 
SCRIPTS  = $(wildcard $(SCRIPTDIR)/*.cpp)
BIN      = $(SCRIPTS:$(SCRIPTDIR)/%.cpp=$(BINDIR)/%.x) 
BINDEBUG = $(SCRIPTS:$(SCRIPTDIR)/%.cpp=$(BINDIR)/%_debug.x) 

SRC      = $(wildcard $(SRCDIR)/*.cpp)
OBJ      = $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o) 
OBJDEBUG = $(SRC:$(SRCDIR)/%.cpp=$(OBJDIR)/%_debug.o)
	
release: $(BIN)
debug:   $(BINDEBUG)

all: release debug

############# Release version

# source files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(IPATH)

#script source files
$(SCRIPTOBJDIR)/%.o: $(SCRIPTDIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(IPATH)

# Executables
$(BINDIR)/%.x: $(SCRIPTOBJDIR)/%.o $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(IPATH) $(LPATH) $(LIBS)

############## Debugging version

# source files
$(OBJDIR)/%_debug.o: $(SRCDIR)/%.cpp
	$(CXXDEBUG) $(CXXFLAGS) -c $< -o $@ $(IPATH)

#script source files
$(SCRIPTOBJDIR)/%_debug.o: $(SCRIPTDIR)/%.cpp
	$(CXXDEBUG) $(CXXFLAGS) -c $< -o $@ $(IPATH)

# Executables
$(BINDIR)/%_debug.x: $(SCRIPTOBJDIR)/%_debug.o $(OBJ)
	$(CXXDEBUG) $(CXXFLAGS) -o $@ $^ $(IPATH) $(LPATH) $(LIBS)

.PHONY: clean
clean:
	-rm -f $(OBJDIR)/*.o;
	-rm -f $(SCRIPTOBJDIR)/*.o;
	-rm -f $(BINDIR)/*.x;
