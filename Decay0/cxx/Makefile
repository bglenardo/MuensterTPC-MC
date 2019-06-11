ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs) -lm -ldl
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs) -lm -ldl
INCLUDES      = $(shell $(ROOTSYS)/bin/root-config --cflags)

#ARCH          = -m32
OPT           = -g -march=core2 -mfpmath=sse -fopenmp
CXX           = g++
CXXFLAGS      = -Wall $(OPT) $(ARCH) -fPIC
INCS          =  $(INCLUDES) 
LD            = g++
LDFLAGS       = $(OPT)
SOFLAGS       = -Wl,-soname -shared
LIB_EXTRA     =
LIBS          = $(ROOTGLIBS) $(LIB_EXTRA) -L/usr/X11R6/lib -lXpm -lX11 
%.o:    %.cxx
	$(CXX) -c $(INCS) $(CXXFLAGS) $<

%.o:    %.C
	$(CXX) -c $(INCS) $(CXXFLAGS) $<
Decay:   Decay.o
	$(LD) $(LDFLAGS) -o Decay Decay.o $(LIBS) 

clean: 
	rm Decay.o
