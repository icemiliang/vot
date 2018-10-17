# vot: Variational Optimal Transportation
#
# Author   : Mi,Liang (Arizona State University)
# Email    : icemiliang@gmail.com
# Date     : Sept 28th 2018

# C++ compiler
CXX=g++

# Diractories
LIBDIR_DEMO = ../lib/
INCDIR_DEMO = ../include/
SRCDIR_DEMO = ../src/vot/
SRCDIRX_DEMO = ../src/votx/

LIBDIR = ../../lib/
INCDIR = ../../include/
SRCDIR = ../../src/vot/
SRCDIRX = ../../src/votx/

# Flags for the C++ compiler
CFLAGS=-Wall -ansi -pedantic -O3 -std=c++11

# Include and library paths
E_INC=-I$(SRCDIR) -I$(SRCDIR_DEMO) -I$(SRCDIRX) -I$(SRCDIRX_DEMO) -I$(INCDIR) -I$(INCDIR_DEMO) -I/usr/include/
E_LIB=-L$(LIBDIR) -L$(LIBDIR_DEMO) -L/usr/lib/
