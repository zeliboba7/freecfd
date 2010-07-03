src = variable.cc

.SUFFIXES: .cc .o

objects =  $(src:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

inputs : $(objects)
	ar rcs libvariable.a $(objects)

.PHONY : clean install

clean: 
	-rm $(objects)

install:
	cp variable.h ../../include
	mv libvariable.a ../../lib
