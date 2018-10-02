# vot makefile
#
# Author   : Mi,Liang (Arizona State University)
# Email    : icemiliang@gmail.com
# Date     : Sept 27th 2018

# Load the configuration file
include ../config.mk

# List of the source files
objs=cell.o common.o container.o unitcell.o v_compute.o c_loops.o \
     v_base.o wall.o container_prd.o diagram.o ot.o
src=$(patsubst %.o,%.cc,$(objs))

# Makefile rules
all: $(objs) libvot.a

depend:
	$(CXX) $(E_INC) -MM $(src) >Makefile.dep

include Makefile.dep

libvot.a: $(objs)
	rm -f libvot.a
	ar rs libvot.a $^
	mv libvot.a $(LIBDIR)

%.o: %.cc
	$(CXX) $(CFLAGS) $(E_INC) -c $<

help: Doxyfile $(SOURCE)
	doxygen Doxyfile

clean:
	rm -f $(objs)
	rm -f $(LIBDIR)/libvot.a
	