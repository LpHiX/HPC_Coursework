.PHONY : default run clean doc debug

CXX = g++
CXXFLAGS = -Werror -Wall -O2 -ftree-vectorize -I$(SRCDIR)
LIBS = -llapack -lblas -lboost_program_options
SRCDIR = ./src
TARGET = poisson

SOURCES = $(wildcard $(SRCDIR)/*.cpp) 
OBJS = $(SOURCES:.cpp=.o)
HDRS = $(wildcard $(SRCDIR)/*.h)

default: run

debug: CXXFLAGS = -g -O0 -Wall -Werror -I$(SRCDIR)
debug: clean $(TARGET)

# Compiles all .cpp to .o
$(SRCDIR)/%.o : $(SRCDIR)/%.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ -c $< 

# Links
$(TARGET) : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)
	rm -f $(SRCDIR)/*.o 
	
run: $(TARGET)
	@echo "\033[0;32mRunning Poisson \033[0m"
	./$(TARGET)

clean:
	rm -f $(SRCDIR)/*.o 

doc:
	doxygen Doxyfile
	@echo "\033[0;34mDocumentation generated in /docs folder.\033[0m"