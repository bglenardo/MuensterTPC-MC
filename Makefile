# --------------------------------------------------------------
# GNUmakefile of MuensterTPCsim
#
# @author Lutz Althueser
# @date   2015-07-07
#
# changes:
#		20150928: - added support for changeable main class file name
#							- added 'link' parameter to automatic link the binaries
#							- disable warnings with the '-w' flag
#							- turn on optimmization by deleting the '-g' flag
# --------------------------------------------------------------

# get full name of the main class file

name := MuensterTPC-MC
G4TARGET := $(name)
G4EXLIB := true
MYEXEPATH:=$(G4WORKDIR)/bin/Darwin-clang
#CPPVERBOSE := true

G4DEBUG := 0

GEANTLIBS       = $(shell geant4-config --libs)
ROOTCFLAGS      = $(shell root-config --cflags) -Wno-shadow -w
ROOTLIBS        = $(shell root-config --nonew --libs)
ROOTGLIBS       = $(shell root-config --glibs)

EXTRALIBS +=$(ROOTLIBS) $(GEANTLIBS)
CPPFLAGS += $(ROOTCFLAGS)

.PHONY: all
all: lib bin

#change because of compiler error
include $(G4INSTALL)/config/alt_binmake.gmk
#include $(G4INSTALL)/config/binmake.gmk


# call this routine with 'make link' to create a new symlink of the binary
link: 
	[ -f $(name) ] || ln -s $(G4WORKDIR)/bin/$(G4SYSTEM)/$(name) ./$(name)
