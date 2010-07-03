src = vec3d.cc

.SUFFIXES: .cc .o

objects =  $(src:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

vec3d : $(objects)
	ar rcs libvec3d.a $(objects)

.PHONY : clean install

clean: 
	-rm $(objects)

install:
	cp vec3d.h ../../include
	mv libvec3d.a ../../lib
