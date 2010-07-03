src = material.cc

.SUFFIXES: .cc .o

objects =  $(src:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

material : $(objects)
	ar rcs libmaterial.a $(objects)

.PHONY : clean install

clean: 
	-rm $(objects)

install:
	cp material.h ../../include
	mv libmaterial.a ../../lib
