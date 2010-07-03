src = polynomial.cc

.SUFFIXES: .cc .o

objects =  $(src:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

polynomial : $(objects)
	ar rcs libpolynomial.a $(objects)

.PHONY : clean install

clean: 
	-rm $(objects)

install:
	cp polynomial.h ../../include
	mv libpolynomial.a ../../lib
