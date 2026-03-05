.PHONY : default run clean doc debug test

INCDIR = include
SRCDIR = src
BUILDDIR = build
TESTDIR = tests


CXX = g++
CXXFLAGS = -O2 -Werror -Wall -ftree-vectorize -I $(INCDIR)
# LIBS = -llapack -lblas -lboost_program_options 
LIBS = -lboost_timer -lboost_chrono -lboost_system -lboost_program_options -llapack -lblas -lrt
TARGET = poisson
TEST_TARGET = run_tests



SRCS = $(wildcard $(SRCDIR)/*.cpp) 
TEST_SRCS = $(wildcard $(TESTDIR)/*.cpp) 

# HDRS = $(wildcard $(INCDIR)/*.h)

OBJS = $(patsubst $(SRCDIR)/%.cpp, $(BUILDDIR)/%.o, $(SRCS))
TEST_OBJS = $(patsubst $(TESTDIR)/%.cpp, $(BUILDDIR)/%.o, $(TEST_SRCS))

MAIN_OBJ = $(BUILDDIR)/$(TARGET).o
SHARED_OBJS = $(filter-out $(MAIN_OBJ), $(OBJS))

default: run
# default: run test

debug: CXXFLAGS = -g -O0 -Wall -Werror -I$(INCDIR)
test: CXXFLAGS =  -O0 -Wall -Werror -ftree-vectorize -I$(INCDIR)
debug: clean $(TARGET)

$(BUILDDIR):
	mkdir -p $@

# Compiles all .cpp to .o
$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $< 

$(BUILDDIR)/%.o : $(TESTDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $< 

# Links
$(TARGET) : $(MAIN_OBJ) $(SHARED_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

$(TEST_TARGET) : $(TEST_OBJS) $(SHARED_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

run: $(TARGET)
	@echo "\033[0;32mRunning Poisson \033[0m"
	./$(TARGET)

test: $(TEST_TARGET)
	@echo "\033[0;32mRunning Tests \033[0m"
	./$(TEST_TARGET)

clean:
	rm -rf $(BUILDDIR) $(TARGET) $(TEST_TARGET)

doc:
	doxygen Doxyfile
	@echo "\033[0;34mDocumentation generated in /docs folder.\033[0m"