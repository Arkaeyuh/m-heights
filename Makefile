CXX       := g++
CXXFLAGS  := -std=c++17 -O3 -march=native -fopenmp -Wall -Wextra -pedantic
INCLUDES  := -Iinclude
LDFLAGS   := -lglpk -fopenmp

SRCDIR    := src
BUILDDIR  := build
TARGET    := combined_annealing

SRC       := $(wildcard $(SRCDIR)/*.cpp)
OBJ       := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(SRC))

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

clean:
	rm -rf $(BUILDDIR) $(TARGET)
