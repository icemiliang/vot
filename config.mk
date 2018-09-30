# vot: Variational Optimal Transportation
#
# Author   : Mi,Liang (Arizona State University)
# Email    : icemiliang@gmail.com
# Date     : Sept 28th 2018

# C++ compiler
CXX=g++

# Diractories
LIBDIR = ../lib/
INCDIR = ../include/
SRCDIR = ../src/

# Flags for the C++ compiler
CFLAGS=-Wall -ansi -pedantic -O3 -std=c++11

# Include and library paths for compilation of the examples
E_INC=-I$(SRCDIR) -I$(INCDIR) -I/usr/include/
E_LIB=-L$(SRCDIR) -L$(LIBDIR) -L/usr/lib/
