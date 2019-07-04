.PHONY: all clean

.DEFAULT: all


CXX = c++

CXXFLAGS = -ansi -pedantic -std=c++11 -fopenmp

TARGETS = run_quantumsim

all: $(TARGETS)

clean:
	rm -f $(TARGETS)

$(TARGETS) : %: main_quantumsim.cc
	$(CXX) $(CXXFLAGS) -o $@ $<
