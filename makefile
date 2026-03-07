.PHONY : default run clean doc debug test

INCDIR = include
SRCDIR = src
BUILDDIR = build
TESTDIR = tests


CXX = g++
poisson-mpi : CXX = mpicxx
poisson-cuda : CXX = nvcc

CXXFLAGS = -O3 -Werror -Wall -ftree-vectorize  -I $(INCDIR)
# LIBS = -llapack -lblas -lboost_program_options 
LIBS = -lboost_timer -lboost_chrono -lboost_system -lboost_program_options -llapack -lblas -lrt
TARGETS = poisson poisson-mpi 
# TARGETS = poisson poisson-mpi poisson-cuda
TEST_TARGET = run_tests



SRCS = $(wildcard $(SRCDIR)/*.cpp) 
TEST_SRCS = $(wildcard $(TESTDIR)/*.cpp) 

# HDRS = $(wildcard $(INCDIR)/*.h)

OBJS = $(patsubst $(SRCDIR)/%.cpp, $(BUILDDIR)/%.o, $(SRCS))
TEST_OBJS = $(patsubst $(TESTDIR)/%.cpp, $(BUILDDIR)/%.o, $(TEST_SRCS))

ALL_MAINS = $(patsubst %, $(BUILDDIR)/%.o, $(TARGETS))
SHARED_OBJS = $(filter-out $(ALL_MAINS), $(OBJS))

default: runmpi
# default: run test

# debug: CXXFLAGS = -g -O0 -Wall -Werror -I$(INCDIR)
# test: CXXFLAGS =  -O0 -Wall -Werror -ftree-vectorize -I$(INCDIR)
# debug: clean $(TARGETS)

debug:
	$(MAKE) clean
	$(MAKE) $(TARGETS) CXXFLAGS="-g -O0 -Wall -Werror -I$(INCDIR)"

optimize:
	$(MAKE) clean
	$(MAKE) $(TARGETS) CXXFLAGS="-g -O3 -ftree-vectorize -Wall -Werror -I$(INCDIR)"

$(BUILDDIR):
	mkdir -p $@

# Compiles all .cpp to .o

$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $< 

$(BUILDDIR)/%.o : $(SRCDIR)/%.cu | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(BUILDDIR)/%.o : $(TESTDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $< 

# Links
$(TARGETS) : % :  $(BUILDDIR)/%.o $(SHARED_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(TEST_TARGET) : $(TEST_OBJS) $(SHARED_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

run: poisson
	@echo "\033[0;32mRunning Poisson \033[0m"
	./poisson

NP = 12
runmpi: poisson-mpi
	@echo "\033[0;32mRunning Poisson \033[0m"
	mpiexec --mca btl_tcp_if_include lo -np $(NP) ./poisson-mpi --Px 1 --Py 1 --Pz 1

tests: $(TEST_TARGET)
	@echo "\033[0;32mRunning Tests \033[0m"
	./$(TEST_TARGET)



clean:
	rm -rf $(BUILDDIR) $(TARGET) $(TEST_TARGET)

doc:
	doxygen Doxyfile
	@echo "\033[0;34mDocumentation generated in /docs folder.\033[0m"

cleandoc :
	rm -rf docs