# vot makefile
#
# Author   : Mi,Liang (Arizona State University)
# Email    : icemiliang@gmail.com
# Date     : Oct 16th 2018

# Load the configuration file
include ../../config.mk

# List of the source files
objs=otx.o
src=$(patsubst %.o,%.cc,$(objs))

# Makefile rules
all: $(objs) libvotx.a 

depend:
	$(CXX) $(E_INCX) -MM $(src) >Makefile.dep

include Makefile.dep

libvotx.a: $(objs)
	rm -f libvotx.a
	ar rs libvotx.a $^
	mv libvotx.a $(LIBDIR)

%.o: %.cc
	$(CXX) $(CFLAGS) $(E_INCX) -c $<

clean:
	rm -f $(objs)
	rm -f $(LIBDIR)/libvotx.a


