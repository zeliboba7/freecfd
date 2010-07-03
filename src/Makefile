src = \
bc_interface.cc \
bc_interface_sync.cc \
main.cc \
read_inputs.cc \
read_restart.cc \
set_bcs.cc \
time_step.cc \
write_output.cc \
write_restart.cc 

.SUFFIXES: .cc .o

objects =  $(src:.cc=.o)

.cc.o:
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $<

fcfd : $(objects)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS) -o $@ $(objects)

.PHONY : clean install

clean: 
	-rm $(objects) fcfd

install:
	mv fcfd ../bin


