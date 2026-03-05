.PHONY : default run clean doc debug test

INCDIR = include
SRCDIR = src
BUILDDIR = build
TESTDIR = tests


CXX = g++
CXXFLAGS = -Werror -Wall -ftree-vectorize -I $(INCDIR)
LIBS = -llapack -lblas -lboost_program_options
TARGET = poisson
TEST_TARGET = run_tests



SRCS = $(wildcard $(SRCDIR)/*.cpp) 
TEST_SRCS = $(wildcard $(TESTDIR)/*.cpp) 

# HDRS = $(wildcard $(INCDIR)/*.h)

OBJS = $(patsubst $(SRCDIR)/%.cpp, $(BUILDDIR)/%.o, $(SRCS))
TEST_OBJS = $(patsubst $(TESTDIR)/%.cpp, $(BUILDDIR)/%.o, $(TEST_SRCS))

MAIN_OBJ = $(BUILDDIR)/$(TARGET).o
SHARED_OBJS = $(filter-out $(MAIN_OBJ), $(OBJS))

default: run test

debug: CXXFLAGS = -g -O0 -Wall -Werror -I$(INCDIR)
debug: clean $(TARGET)

$(BUILDDIR):
	mkdir -p $@

# Compiles all .cpp to .o
$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -O2 -o $@ -c $< 

$(BUILDDIR)/%.o : $(TESTDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -O0 -o $@ -c $< 

# Links
$(TARGET) : $(MAIN_OBJ) $(SHARED_OBJS)
	$(CXX) $(CXXFLAGS) -O2 -o $@ $^ $(LIBS)

$(TEST_TARGET) : $(TEST_OBJS) $(SHARED_OBJS)
	$(CXX) $(CXXFLAGS) -O0 -o $@ $^ $(LIBS)

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