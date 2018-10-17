# vot makefile
#
# Author   : Mi,Liang (Arizona State University)
# Email    : icemiliang@gmail.com
# Date     : Sept 27th 2018

# Load the configuration file
include ../config.mk

# List of executables
EXECUTABLES= vot votx

# Makefile rules
all: $(EXECUTABLES)

vot: vot.cpp
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o vot vot.cpp -lvot -lboost_filesystem -lboost_program_options -lboost_system

votx: votx.cpp
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o votx votx.cpp -lvotx -lboost_filesystem -lboost_program_options -lboost_system

clean:
	rm -f $(EXECUTABLES)

