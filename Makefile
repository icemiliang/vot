# vot makefile
#
# Author   : Mi,Liang (Arizona State University)
# Email    : icemiliang@gmail.com
# Date     : Sept 27th 2018

# Tell make that these are phony targets
.PHONY: all help clean

include config.mk

# Build all of the executable files
all:
	$(MAKE) -C src
	$(MAKE) -C demo

# Build the help files (with Doxygen)
help:
	$(MAKE) -C src help

# Clean up the executable files
clean:
	$(MAKE) -C src clean
	$(MAKE) -C demo clean
