# vot makefile
#
# Author   : Mi,Liang (Arizona State University)
# Email    : icemiliang@gmail.com
# Date     : Oct 16th 2018

include config.mk

all:
	$(MAKE) -C src
	$(MAKE) -C demo

clean:
	$(MAKE) -C src clean
	$(MAKE) -C demo clean
