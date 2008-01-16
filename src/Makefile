CC=mpic++ 
#OPTS=-Wall -g
OPTS=-O3
EXE=fcfd

INCLUDES= -I$(CGNS_BASE)/include -I$(PARMETIS_BASE)/include
LIBS= -L$(CGNS_BASE)/lib -L$(PARMETIS_BASE)/lib

SRC= \
	main.cc grid.cc vec3d.cc RoeFlux.cc  write_tec.cc writeTecplotMacro.cc inputs.cc sparse.cc readInputs.cc initialize.cc setBCs.cc fou.cc read_tec.cc hancock.cc diffFlux.cc
.SUFFIXES: .cc .o

OBJ=	$(SRC:.cc=.o)

.cc.o:
	$(CC) $(OPTS) $(INCLUDES) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(OPTS) -o $@ $(OBJ) $(LIBS) -lcgns -lparmetis -lmetis

clean:	
	rm -rf $(OBJ) $(EXE)
