# You can rebuild fortran sources from their original git archive as
# follows
#
#   git clone https://github.com/cepek/geodet.pc.git
#
# and run make
#
#
 all : mtgso1.f mtgso2.f

mtgso1.f : geodet.pc/geolib.pc/mt/mtgso1.for Makefile gso2022.sed
	sed -f gso2022.sed < $< >  mtgso1.f

mtgso2.f : geodet.pc/geolib.pc/mt/mtgso2.for Makefile gso2022.sed
	sed -f gso2022.sed < $< >  mtgso2.f

clean:
	rm -f *~
