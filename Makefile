# A Generic Makefile for compiling ROOT programs
# R. Michaels, rom@jlab.org, Aug 2001  See also README !!
# Version of this release
# 
ROOTLIBS      = $(shell root-config --libs) -lMathMore
ROOTGLIBS     = $(shell root-config --glibs)
LHAPDFLIBS    = 
CXX           = g++
CXXFLAGS      = -lgfortran -lgsl -lgslcblas -lm -Wall -frtti -fexceptions -fPIC \
                   -DLINUXVERS -I$(ROOTSYS)/include -O

FFLAGS	      = -lgfortran -O -W -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8

# Linux with g++
INCLUDES      = -I$(ROOTSYS)/include 
CXX           = g++ -fopenmp
LD            = g++
LDFLAGS       = -O2

LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)

ALL_LIBS =  $(GLIBS) $(LIBS) $(LHAPDFLIBS)

# The following sources comprise the package of scaler classes by R. Michaels.
# O = main
SRC = $(O).C

HEAD = $(SRC:.C=.h)
DEPS = $(SRC:.C=.d)
SCALER_OBJS = $(SRC:.C=.o)

F_OBJS = 
cc_OBJS =
total_OBJS = $(F_OBJS) $(cc_OBJS)

# Test code executibles
PROGS = $(O)

%.o: %.f
	$(CXX) $(FFLAGS) -c $< -o $@

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(O): $(O).o $(O).C $(total_OBJS)
	rm -f $@
	$(CXX) $(CXXFLAGS) -o $@ $(O).o $(total_OBJS) $(ALL_LIBS)


clean:
	rm -f *.o core *~ *.d *.tar $(PROGS)

realclean:  clean
	rm -f *.d

###

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .C .o .d

%.o:	%.C
	$(CXX) $(CXXFLAGS) -c $<

%.d:	%.C
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) -c $< \
		| sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

-include $(DEPS)








