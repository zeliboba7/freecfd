src = utilities.cc

.SUFFIXES: .cc .o

objects =  $(src:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

utilities : $(objects)
	ar rcs libutilities.a $(objects)

.PHONY : clean install

clean: 
	-rm $(objects)

install:
	cp utilities.h ../../include
	mv libutilities.a ../../lib
