.PHONY : default run clean

CXX = g++
CXXFLAGS = -llapack -lblas -Werror

default: run clean

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<

poisson : poisson.cpp
	$(CXX) $(CXXFLAGS) -o $@ $<
	

run: poisson
	@echo "\033[0;32mRunning Poisson \033[0m"
	./poisson

clean:
	rm -f poisson
	