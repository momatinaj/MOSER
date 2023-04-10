ESCAPE_HOME := .

OBJECTS := Graph.o GraphIO.o TriangleProgram.o

TARGETS := libescape.a

all: exe

libescape.a : $(OBJECTS)
	ar cruv $@ $^

include common.mk


exe: libescape.a
	$(MAKE) -C exe

cleanexe:
	$(MAKE) -C exe clean cleandep

.PHONY: exe



