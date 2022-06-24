# sed '/^#/d'
#
INPUTS  := $(wildcard *.dat)
OUTPUTS := $(patsubst %.dat,%.txt,$(INPUTS))

all : $(OUTPUTS)

%.txt : %.dat
	sed '/^#/d' $<  > $@
