CC=mpic++ 
#OPTS=-Wall -g
OPTS=-O3
EXE=fcfd

INCLUDES= /home/esozer/cgnslib_2.5/include/ -I/usr/include/parmetis/ 
LIBS= /home/esozer/cgnslib_2.5/lib/libcgns.a -L/usr/lib  -lparmetis -lmetis  -lm

SRC= \
	main.cc grid.cc vec3d.cc RoeFlux.cc  write_tec.cc writeTecplotMacro.cc inputs.cc sparse.cc readInputs.cc initialize.cc setBCs.cc fou.cc read_tec.cc hancock.cc diffFlux.cc
.SUFFIXES: .cc .o

OBJ=	$(SRC:.cc=.o)

.cc.o:
	$(CC) $(OPTS) -I$(INCLUDES) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(OPTS) -o $@ $(OBJ) $(LIBS)

clean:	
	rm -rf $(OBJ) $(EXE)
