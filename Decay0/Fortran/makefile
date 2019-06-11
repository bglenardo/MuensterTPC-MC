#Lutz Althüser
#Compilieren
#FCC=xlf
FCC=gfortran

#Compiler-Flags
#FFLAGS += -qextname -g
#FFLAGS += -g
FFLAGS += -m32
#Linken
LD=gfortran

#Linker-Flags
LDFLAGS= -m32

#Libraries
#LIBS += -L/usr/local/cern/2006/lib \
# -lmathlib -lkernlib #-lpacklib

LIBS += -lmathlib -lkernlib -lpacklib
#Objektfiles
#OBJS = espectrum.o
#OBJS = .o
#OBJS = sensi.o
#OBJS = atome.o
####OBJS = decay0.o
OBJS = decay0.o
#OBJS = hwzcheck.o
#OBJS = hwzchecknew.o
#OBJS = naglib.o
#OBJS = bbspec.o

# Name des Programms
#PROGRAM = espectrum
#PROGRAM = sensi
PROGRAM = decay0
#PROGRAM = decay0.1.new
#PROGRAM = atome
#PROGRAM = hwzcheck
#PROGRAM = hwzchecknew
#PROGRAM = bbspec

# siehe naechster Punkt
all: $(PROGRAM)

# Linken
$(PROGRAM): $(OBJS)
	$(LD) $(OBJS) $(LIBS) -o $(PROGRAM)


# Compilieren
%.o: %.for
	$(FCC) -c -o $@ $<


.PHONY :
clean:
	$(RM) *.o
