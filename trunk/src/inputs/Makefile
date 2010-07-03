src = inputs.cc

.SUFFIXES: .cc .o

objects =  $(src:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

inputs : $(objects)
	ar rcs libinputs.a $(objects)

.PHONY : clean install

clean: 
	-rm $(objects)

install:
	cp inputs.h ../../include
	mv libinputs.a ../../lib
