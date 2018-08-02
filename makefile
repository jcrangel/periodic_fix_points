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
CFLAGS := -w -O2 # -Wall
INC := -I/usr/local/boost_1_67_0/ -I/usr/local/eigen 

all: $(TARGET)

$(TARGET): FORCE
	@echo " CompLinking..."
	@echo " $(CC) $(CFLAGS) $(SOURCES) -o $(TARGET) $(INC)";  $(CC) $(CFLAGS) $(SOURCES) -o $(TARGET) $(INC)

FORCE:

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)
	
.PHONY: clean