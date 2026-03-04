.PHONY : default run clean docs

CXX = g++
CXXFLAGS = -Werror -Wall -O2 -ftree-vectorize
LIBS = -llapack -lblas
srcdir = ./src

default: run clean

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

poisson : $(srcdir)/poisson.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
	

run: poisson
	@echo "\033[0;32mRunning Poisson \033[0m"
	./poisson

clean:
	rm -f poisson
	
docs:
	doxygen Doxyfile
	@echo "\033[0;34mDocumentation generated in /docs folder.\033[0m"