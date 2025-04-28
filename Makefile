# Makefile for combined_annealing

# C++ compiler and flags
CXX       := g++
CXXFLAGS  := -std=c++17 -O3 -march=native -fopenmp -Wall -Wextra -pedantic
INCLUDES  := -Iinclude

# Linker flags
LDFLAGS   := -lglpk

# Directories
SRCDIR    := src
BUILDDIR  := build
TARGET    := combined_annealing

# Source files
SRC       := $(wildcard $(SRCDIR)/*.cpp)
OBJ       := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(SRC))

.PHONY: all clean

all: $(TARGET)

# Link step
$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $^ $(LDFLAGS)

# Compile step
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Run the program 
run: $(TARGET)
	./$(TARGET)

# Clean up
clean:
	rm -rf $(BUILDDIR) $(TARGET)
