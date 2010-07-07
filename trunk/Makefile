export CXX = mpic++
export CXXFLAGS = -O3

FCFD_HOME = /Users/esozer/fcfd
CGNS_HOME = /Users/esozer/fcfd/deps/cgnslib_2.5
PARMETIS_HOME = /Users/esozer/fcfd/deps/ParMetis-3.1.1
PETSC_HOME=/Users/esozer/fcfd/deps/petsc-3.1p2

export INCLUDE = -I$(FCFD_HOME)/src -I$(FCFD_HOME)/include -I$(CGNS_HOME)/include -I$(PARMETIS_HOME)/include -I$(PETSC_HOME)/include
export LIBS = -L$(FCFD_HOME)/lib -L$(CGNS_HOME)/lib -L$(PARMETIS_HOME)/lib -L$(PETSC_HOME)/lib -lutilities -lvec3d -lpolynomial -linputs -lmaterial -lgrid -lvariable -lns -lhc -lcgns -lmetis -lparmetis -lpetsc

SUBDIRS= src/vec3d src/polynomial src/utilities src/inputs src/grid src/variable src/material src/navier_stokes src/heat_conduction src

.PHONY: subdirs $(SUBDIRS) clean

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@
	$(MAKE) -C $@ install

clean: 
	-$(RM) -f ./bin/*
	-$(RM) -f ./lib/*
	-$(RM) -f ./include/*
	-for d in $(SUBDIRS); do (cd $$d; $(MAKE) clean ); done
