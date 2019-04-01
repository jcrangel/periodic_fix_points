#
# Simple make file for one main.cpp with multiples .h
#FORCE , force to always make the project
#
 
CC := g++ # This is the main compiler
SRCDIR := src
BUILDDIR := build
TARGET := bin/runner
 
SRCEXT := cpp
SOURCES := src/main.cpp
CFLAGS := -w -O2 -std=c++11# -Wall
INC := -I/home/rangel/boost_1_69_0/ -I/home/rangel/eigen 

all: $(TARGET)

$(TARGET): FORCE
	@echo " CompLinking..."
	@echo " $(CC) $(CFLAGS) $(SOURCES) -o $(TARGET) $(INC)";  $(CC) $(CFLAGS) $(SOURCES) -o $(TARGET) $(INC)

FORCE:

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)
	
.PHONY: clean