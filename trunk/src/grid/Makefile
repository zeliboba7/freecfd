src = \
grid.cc \
grid_create_elements.cc \
grid_readers.cc \
grid_averaging_metrics.cc \
grid_partition.cc \
grid_transform.cc 

.SUFFIXES: .cc .o

objects =  $(src:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

grid : $(objects)
	ar rcs libgrid.a $(objects)

.PHONY : clean install

clean: 
	-rm $(objects)

install:
	cp grid.h ../../include
	mv libgrid.a ../../lib
